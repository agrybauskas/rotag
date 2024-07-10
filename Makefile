#
# C++ and C libraries and programs.
#

BIN_DIR=bin
SRC_DIR=src

LIB_DIR=${SRC_DIR}/lib
LIB_SRC=${wildcard ${LIB_DIR}/*.cpp ${LIB_DIR}/ForceField/*.cpp}
OBJ_DIR=${SRC_DIR}/lib

BIN_SRC=$(wildcard ${SRC_DIR}/*.cpp)
HEADERS=${LIB_SRC:%.cpp=%.h}

PARSER=$(wildcard ${LIB_DIR}/Grammar/*.y)
PARSER_SRC=$(PARSER:%.y=%.cpp)
PARSER_HEADERS=${PARSER_SRC:%.cpp=%.h}

LEXER=$(wildcard ${LIB_DIR}/Grammar/*.l)
LEXER_SRC=$(LEXER:%.l=%.cpp)
LEXER_HEADERS=${LEXER_SRC:%.cpp=%.h}

CPP_OBJS=${LIB_SRC:%.cpp=%.o}
CPP_OBJS+=${PARSER_SRC:%.cpp=%.o}
CPP_OBJS+=${LEXER_SRC:%.cpp=%.o}

CPP_BIN=${BIN_SRC:${SRC_DIR}/%.cpp=${BIN_DIR}/%}
CPP_LIB=-lboost_filesystem

C_LIBDIR=-Isrc/externals/cexceptions -Isrc/externals/codcif -Isrc/externals/getoptions
C_OBJS=${SRC_DIR}/externals/codcif/obj/*.o ${SRC_DIR}/externals/cexceptions/obj/*.o ${SRC_DIR}/externals/getoptions/obj/*.o

.PHONY: all
.PRECIOUS: ${CPP_OBJS} ${PARSER_SRC} ${LEXER_SRC} ${PARSER_HEADERS} ${LEXER_HEADERS}

all: build-externals | ${CPP_BIN}

build-externals:
	make -C src/externals/cexceptions
	make -C src/externals/codcif
	make -C src/externals/getoptions

%.cpp: %.l
	flex --header-file=$(basename $@).h -o $@ $<

%.cpp: %.y
	bison --defines=$(basename $@).h -o $@ $<

%.o: %.cpp %.h
	g++ -c -Wall -std=c++11 -g -o $@ $< ${CPP_LIB} ${C_LIBDIR}

${BIN_DIR}/%: ${SRC_DIR}/%.cpp ${CPP_OBJS}
	g++ -Wall -std=c++11 -g -o $@ $< ${CPP_OBJS} ${C_OBJS} ${CPP_LIB} ${C_LIBDIR}

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
	rm -f ${CPP_OBJS}
	rm -f ${CPP_BIN}
	rm -f ${LIB_DIR}/Grammar/*.cpp
	rm -f ${LIB_DIR}/Grammar/*.h
	rm -f ${LIB_DIR}/Grammar/*.o
