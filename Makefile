all:

#
# Grammar compilation.
#

YAPP_DIR=lib/Grammar
YAPP_FILES=${wildcard ${YAPP_DIR}/*.yp}
GRAMMAR_MODULES=${YAPP_FILES:%.yp=%.pm}

${YAPP_DIR}/%.pm: ${YAPP_DIR}/%.yp
	yapp -o $@ $<

#
# Compiling CPP and linking to Perl5 with SWIG.
#

# ${CPP_DIR}/%.o: ${CPP_DIR}/%.cpp
# 	g++ -c -I${CPP_DIR} $< -o $@
# CPP_OBJS=${CPP_FILES:%.cpp=%.o}

LIB_DIR=lib
CPP_DIR=${LIB_DIR}/CPP
SWIG_FILES=${wildcard ${CPP_DIR}/*.i}
CPP_FILES=${SWIG_FILES:%.i=%.cpp}
CPP_OBJS=${SWIG_FILES:%.i=%.o}
CPP_LIBS=-lboost_regex
PM_FILES=${SWIG_FILES:%.i=%.pm}
WRAP_FILES=${SWIG_FILES:%.i=%_wrap.cxx}
WRAP_OBJS=${WRAP_FILES:%.cxx=%.o}
SHARED_OBJS=${CPP_OBJS:%.o=%.so}

# CPP_TEST_SRC=tests/src
# CPP_TEST_BIN=tests/bin
# CPP_TEST_FILES=${wildcard ${CPP_TEST_SRC}/*.cpp}
# CPP_TEST_BINS=${CPP_TEST_FILES:${CPP_TEST_SRC}/%.cpp=${CPP_TEST_BIN}/%}

.PRECIOUS: ${CPP_OBJS}

%.pm: %.i
	swig -c++ -perl $<

%_wrap.cxx: %.i
	swig -c++ -perl $<

%.o: %.cpp
	g++ -c -fPIC $< -I$$(perl -e 'use Config; print $$Config{archlib};')/CORE -o $@

%_wrap.o: %_wrap.cxx
	g++ -c -fPIC $< -I$$(perl -e 'use Config; print $$Config{archlib};')/CORE \
	    -o $@

%.so: %_wrap.o %.o
	g++ -shared $^ -o $@

.PHONY: all

all: ${GRAMMAR_MODULES} plugins | ${CPP_OBJS} ${PM_FILES} ${WRAP_OBJS} ${SHARED_OBJS}

#
# Instalation of dependencies.
#

.PHONY: build

build:
	./dependencies/$$(lsb_release -is)-$$(lsb_release -rs)/install.sh

#
# Unit tests.
#

TEST_CASES_DIR=tests/cases
TEST_OUT_DIR=tests/outputs
TEST_CASES=${sort ${wildcard ${TEST_CASES_DIR}/*.sh}}
TEST_DIFF=${TEST_CASES:${TEST_CASES_DIR}/%.sh=${TEST_OUT_DIR}/%.diff}

# Common test commands.

define can_run_test
[ ! -e ${TEST_CASES_DIR}/$*.chk ]
endef

.PHONY: test listdiff

test: ${GRAMMAR_MODULES} | ${TEST_DIFF}

listdiff:
	@-find ${TEST_OUT_DIR} -type f -name '*.diff' -size +0 | sort -u

${TEST_OUT_DIR}/%.diff: ${TEST_CASES_DIR}/%.sh ${TEST_OUT_DIR}/%.out
	@if ${can_run_test}; then \
	    ./$< | diff -a -B -w $(basename $@).out - > $@; \
	    if [ $$? -eq 0 ]; \
	    then echo "$<" \
	         | awk '{ printf "%-50s \033[1m[OK]\033[m\n",    $$1 }' \
	         | sed -e 's/ /./g'; \
	    else echo "$<" \
	         | awk '{ printf "%-50s \033[1m[ERROR]\033[m\n", $$1 }' \
	         | sed -e 's/ /./g'; \
                   cat $@; \
	    fi \
	else \
	    echo "$<" \
	        | awk '{ printf "%-50s \033[1m[SKIP]\033[m ", $$1 }' \
	        | sed -e 's/ /./g'; \
	    ${TEST_CASES_DIR}/$*.chk; \
	    touch $@; \
	fi

#
# Plugins
#

PLUGIN_DIR=plugins
PLUGINS=${PLUGIN_DIR}/pymol2-rotag
PLUGINS_ZIP=${PLUGINS:%=%.zip}

.PHONY: plugins

plugins: ${PLUGINS_ZIP}

%.zip: %
	zip -r $@ $<

#
# Coverage.
#

COVERAGE_CASES_DIR=tests/coverage
COVERAGE_CASES=${TEST_CASES:${TEST_CASES_DIR}/%.sh=${COVERAGE_CASES_DIR}/%.sh}
COVERAGE_OUTS=${TEST_CASES:${TEST_CASES_DIR}/%.sh=${COVERAGE_CASES_DIR}/%.out}

.PHONY: coverage

coverage: ${COVERAGE_CASES_DIR}/cover_db/coverage.html

${COVERAGE_CASES_DIR}/cover_db/coverage.html: ${COVERAGE_OUTS}
	cover
	mv cover_db ${COVERAGE_CASES_DIR}

${COVERAGE_CASES_DIR}/%.out: ${COVERAGE_CASES_DIR}/%.sh
	./$< 2>&1 > $@ || true

${COVERAGE_CASES_DIR}/%.sh: ${TEST_CASES_DIR}/%.sh
	cp $^ $@
	sed -i '4i export PERL5OPT=-MDevel::Cover make test' $@

#
# Utilities.
#

.PHONY: clean cleanAll distclean

clean:
	rm -f ${TEST_DIFF}
	rm -f ${COVERAGE_CASES}
	rm -f ${COVERAGE_OUTS}
	rm -fr ${COVERAGE_CASES_DIR}/cover_db

cleanAll distclean: clean
	rm -f ${GRAMMAR_MODULES}
	rm -f ${PERL_MODULE}
	rm -f ${PERL_FORCE_FIELD_MODULE}
	rm -f ${CPP_OBJS}
	rm -f ${CPP_TEST_BINS}
	rm -f ${SHARED_OBJS}
	rm -f ${WRAP_FILES}
	rm -f ${WRAP_OBJS}
	rm -f ${PLUGINS_ZIP}
