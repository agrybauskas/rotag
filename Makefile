#
# C++ and C libraries and programs.
#

BIN_DIR=bin
SRC_DIR=src

LIB_DIR=${SRC_DIR}/lib
LIB_SRC=${wildcard ${LIB_DIR}/*.cc ${LIB_DIR}/ForceField/*.cc ${LIB_DIR}/Grammar/*.cc}
OBJ_DIR=${SRC_DIR}/lib

BIN_SRC=$(wildcard ${SRC_DIR}/*.cc)
HEADERS=${LIB_SRC:%.cc=%.h}

PARSER=$(wildcard ${LIB_DIR}/Grammar/*.y)
PARSER_SRC=$(PARSER:%.y=%.cc)
PARSER_HEADERS=${PARSER_SRC:%.cc=%.h}
PARSER_OBJS=$(PARSER:%.y=%.o)

LEXER=$(wildcard ${LIB_DIR}/Grammar/*.l)
LEXER_SRC=$(LEXER:%.l=%.cc)
LEXER_HEADERS=${LEXER_SRC:%.cc=%.h}
LEXER_OBJS=$(LEXER:%.l=%.o)

CC_OBJS=${LIB_SRC:%.cc=%.o}
CC_OBJS+=${PARSER_SRC:%.cc=%.o}
CC_OBJS+=${LEXER_SRC:%.cc=%.o}

CC_BIN=${BIN_SRC:${SRC_DIR}/%.cc=${BIN_DIR}/%}
CC_LIB=-lboost_filesystem

C_LIBDIR=-Isrc/externals/cexceptions -Isrc/externals/codcif -Isrc/externals/getoptions
C_OBJS=${SRC_DIR}/externals/codcif/obj/*.o ${SRC_DIR}/externals/cexceptions/obj/*.o ${SRC_DIR}/externals/getoptions/obj/*.o

.PHONY: all
.PRECIOUS: ${CC_OBJS} ${PARSER_SRC} ${LEXER_SRC} ${PARSER_HEADERS} ${LEXER_HEADERS}

all: build-externals | ${CC_BIN}

build-externals:
	make -C src/externals/cexceptions
	make -C src/externals/codcif
	make -C src/externals/getoptions

%.cc %.h: %.l
	flex --header-file=$(basename $@).h -o $@ $<

%.cc %.h: %.y
	bison --header=$(basename $@).h -o $@ $<

%.o: %.cc %.h
	g++ -c -Wall -std=c++11 -g -o $@ $< ${CC_LIB} ${C_LIBDIR}

${BIN_DIR}/%: ${SRC_DIR}/%.cc ${CC_OBJS}
	g++ -Wall -std=c++11 -g -o $@ $< ${CC_OBJS} ${C_OBJS} ${CC_LIB} ${C_LIBDIR}

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

test: ${TEST_DIFF} #| ${GRAMMAR_MODULES}

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

#
# Utilities.
#

.PHONY: clean cleanAll distclean

clean: testclean

testclean:
	rm -f ${TEST_DIFF}

cleanAll cleanall distclean: clean
	rm -f ${CC_OBJS}
	rm -f ${CC_BIN}
	rm -f ${LIB_DIR}/Grammar/*.cc
	rm -f ${LIB_DIR}/Grammar/*.h
	rm -f ${LIB_DIR}/Grammar/*.o
