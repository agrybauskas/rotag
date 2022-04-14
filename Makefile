# Environmental variables.

export PERL5LIB:=${PWD}/lib:${PERL5LIB}
export PATH:=${PWD}/scripts:${PATH}

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
# C++ library and program compilation.
#

BIN_DIR=bin
SRC_DIR=src
LIB_DIR=${SRC_DIR}/lib
OBJ_DIR=${SRC_DIR}/lib
LIB_SRC=${wildcard ${LIB_DIR}/*.cpp ${LIB_DIR}/ForceField/*.cpp}
BIN_SRC=$(wildcard ${SRC_DIR}/*.cpp)
HEADERS=${LIB_SRC:%.cpp=%.h}
CPP_OBJS=${LIB_SRC:%.cpp=%.o}
CPP_BIN=${BIN_SRC:${SRC_DIR}/%.cpp=${BIN_DIR}/%}
CPP_LIB=-lboost_regex -lboost_system -lboost_filesystem

.PRECIOUS: ${CPP_OBJS}

%.o: %.cpp %.h
	g++ -c -o $@ $< ${CPP_LIB}

${BIN_DIR}/%: ${SRC_DIR}/%.cpp ${CPP_OBJS}
	g++ -o $@ $< ${CPP_OBJS} ${CPP_LIB}

.PHONY: all

all: ${GRAMMAR_MODULES} | ${CPP_BIN}

#
# Build rule.
#

.PHONY: build

build: all

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

check: test

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

clean: testclean
	rm -f ${COVERAGE_CASES}
	rm -f ${COVERAGE_OUTS}
	rm -fr ${COVERAGE_CASES_DIR}/cover_db

testclean:
	rm -f ${TEST_DIFF}

cleanAll distclean: clean
	rm -f ${GRAMMAR_MODULES}
	rm -f ${CPP_OBJS}
	rm -f ${CPP_BIN}
