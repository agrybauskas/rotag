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
LIB_SRC=${wildcard ${LIB_DIR}/*.cc ${LIB_DIR}/ForceField/*.cc}
BIN_SRC=$(wildcard ${SRC_DIR}/*.cc)
HEADERS=${LIB_SRC:%.cc=%.h}
CC_OBJS=${LIB_SRC:%.cc=%.o}
CC_BIN=${BIN_SRC:${SRC_DIR}/%.cc=${BIN_DIR}/%}
CC_LIB=-lboost_filesystem
C_LIBDIR=-Isrc/externals/cexceptions -Isrc/externals/codcif -Isrc/externals/getoptions
C_OBJS=${SRC_DIR}/externals/codcif/obj/*.o ${SRC_DIR}/externals/cexceptions/obj/*.o ${SRC_DIR}/externals/getoptions/obj/*.o

.PRECIOUS: ${CC_OBJS}

%.o: %.cc %.h
	g++ -c -Wall -g -o $@ $< ${CC_LIB} ${C_LIBDIR}

${BIN_DIR}/%: ${SRC_DIR}/%.cc ${CC_OBJS}
	g++ -Wall -g -o $@ $< ${CC_OBJS} ${C_OBJS} ${CC_LIB}

.PHONY: all

all: build-externals | ${CC_BIN}

build-externals:
	make -C src/externals/cexceptions
	make -C src/externals/codcif
	make -C src/externals/getoptions

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

tags:
	find src/ -name '*.cc' -or -name '*.c' | xargs etags -a -o ${SRC_DIR}/TAGS

clean: testclean
	rm -f ${COVERAGE_CASES}
	rm -f ${COVERAGE_OUTS}
	rm -fr ${COVERAGE_CASES_DIR}/cover_db

testclean:
	rm -f ${TEST_DIFF}

cleanAll distclean: clean
	rm -f ${GRAMMAR_MODULES}
	rm -f ${CC_OBJS}
	rm -f ${CC_BIN}
