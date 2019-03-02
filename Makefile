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
# Generate Perl modules.
#

LIB_DIR=lib
TOOLS_DIR=tools
PERL_TEMPLATE=${LIB_DIR}/Constants.pmin
PERL_MODULE=${LIB_DIR}/Constants.pm

${PERL_MODULE}: ${PERL_TEMPLATE}
	sed 's/@PI@/'$$(${TOOLS_DIR}/calculate-pi)'/g' $^ \
	    | sed 's/@EPSILON@/'$$(${TOOLS_DIR}/calculate-epsilon)'/g' \
	    > $@

#
# Compiling CPP and linking to Perl5 with SWIG.
#

CPP_DIR=${LIB_DIR}
CPP_OBJS=${SWIG_FILES:%.i=%.o}
SWIG_FILES=${wildcard ${CPP_DIR}/*.i}
PM_FILES=${SWIG_FILES:%.i=%.pm}
WRAP_FILES=${SWIG_FILES:%.i=%_wrap.cxx}
WRAP_OBJS=${WRAP_FILES:%.cxx=%.o}
SHARED_OBJS=${CPP_OBJS:%.o=%.so}

CPP_TEST_SRC=tests/src
CPP_TEST_BIN=tests/bin
CPP_TEST_FILES=${wildcard ${CPP_TEST_SRC}/*.cpp}
CPP_TEST_BINS=${CPP_TEST_FILES:${CPP_TEST_SRC}/%.cpp=${CPP_TEST_BIN}/%}

${CPP_DIR}/%.pm: ${CPP_DIR}/%.i
	swig -c++ -perl $<

${CPP_DIR}/%_wrap.cxx: ${CPP_DIR}/%.i
	swig -c++ -perl $<

${CPP_DIR}/%.o: ${CPP_DIR}/%.cpp
	g++ -c -fPIC $< -I$$(perl -e 'use Config; print $$Config{archlib};')/CORE -o $@

${CPP_DIR}/%_wrap.o: ${CPP_DIR}/%_wrap.cxx
	g++ -c -fPIC $< -I$$(perl -e 'use Config; print $$Config{archlib};')/CORE \
	    -o $@

${CPP_DIR}/%.so: ${CPP_DIR}/%_wrap.o ${CPP_DIR}/%.o
	g++ -shared $^ -o $@

${CPP_TEST_BIN}/%: ${CPP_TEST_SRC}/%.cpp
	g++ -c $^ -I$$(perl -e 'use Config; print $$Config{archlib};')/CORE -o $@

#
# Generate force field module.
#

LIB_FORCE_FIELD_DIR=lib/ForceField
TOOLS_DIR=tools
PERL_FORCE_FIELD_TEMPLATE=${LIB_FORCE_FIELD_DIR}/Parameters.pmin
PERL_FORCE_FIELD_CIF=${LIB_FORCE_FIELD_DIR}/parameters.cif
PERL_FORCE_FIELD_MODULE=${LIB_FORCE_FIELD_DIR}/Parameters.pm

${PERL_FORCE_FIELD_MODULE}: ${PERL_FORCE_FIELD_TEMPLATE} ${PERL_FORCE_FIELD_CIF}
	cat $(word 1, $^) > $@
	${TOOLS_DIR}/generate-force-field $(word 2, $^) >> $@
	sed -i 's/\$$/our \$$/g' $@
	sed -i 's/%/our %/g' $@
	sed -i 's/@/our @/g' $@

.PHONY: all

all: ${GRAMMAR_MODULES} ${PERL_MODULE} ${PERL_FORCE_FIELD_MODULE} ${PM_FILES} | ${SHARED_OBJS} ${CPP_TEST_BINS}

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

.PHONY: test

test: ${GRAMMAR_MODULES} ${PERL_MODULE} ${PERL_FORCE_FIELD_MODULE} | ${TEST_DIFF}

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
# Profiler.
#

PROFILER_CASES_DIR=tests/profiler
PROFILER_CASES=${TEST_CASES:${TEST_CASES_DIR}/%.sh=${PROFILER_CASES_DIR}/%.sh}
PROFILER_OUTS=${TEST_CASES:${TEST_CASES_DIR}/%.sh=${PROFILER_CASES_DIR}/%.out}

.PHONY: profiler

profiler: ${PROFILER_CASES_DIR}/nytprof/index.html

${PROFILER_CASES_DIR}/nytprof/index.html: ${PROFILER_OUTS}
	nytprofmerge nytprof.out.*
	nytprofhtml --file nytprof-merged.out
	mv nytprof* ${PROFILER_CASES_DIR}

${PROFILER_CASES_DIR}/%.out: ${PROFILER_CASES_DIR}/%.sh
	./$< 2>&1 > $@ || true

${PROFILER_CASES_DIR}/%.sh: ${TEST_CASES_DIR}/%.sh
	cp $^ $@
	sed -i '4i export PERL5OPT=-d:NYTProf' $@
	sed -i '5i export NYTPROF=addpid=1' $@

#
# Utilities.
#

.PHONY: clean cleanAll distclean

clean:
	rm -f ${TEST_DIFF}
	rm -f ${COVERAGE_CASES}
	rm -f ${COVERAGE_OUTS}
	rm -fr ${COVERAGE_CASES_DIR}/cover_db
	rm -f ${PROFILER_CASES}
	rm -f ${PROFILER_OUTS}
	rm -fr ${PROFILER_CASES_DIR}/nytprof
	rm -fr ${PROFILER_CASES_DIR}/nytprof.out

cleanAll distclean: clean
	rm -f ${GRAMMAR_MODULES}
	rm -f ${PERL_MODULE}
	rm -f ${PERL_FORCE_FIELD_MODULE}
	rm -f ${PM_FILES}
	rm -f ${SHARED_OBJS}
	rm -f ${CPP_OBJS}
	rm -f ${WRAP_FILES}
