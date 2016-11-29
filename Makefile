#
# C++ library build and implementation to Perl5 using SWIG
#

CXX=g++
SRC=lib/cpp

CPP=${wildcard ${SRC}/*.cpp}
SWIG=${wildcard ${SRC}/*.i}

PERL_OBJS=${CPP:%.cpp=%.pm}
SHARED_OBJS=${CPP:%.cpp=%.so}

all: ${SHARED_OBJS} ${PERL_OBJS}

%.pm %_wrap.cxx: %.i
	swig -c++ -perl $<

%.o: %.cpp
	g++ -c $< -o $@

%_wrap.o: %_wrap.cxx
	g++ -c $< `perl -MConfig \
			-e 'print join(" ", \
				       @Config{qw(ccflags optimize cccdlflags)},\
				       "-I$$Config{archlib}/CORE")'` \
	    -o $@

%.so: %.o %_wrap.o
	g++ -shared $^ -o $@

#
# Unit tests
#

TEST_CASES_DIR=test/case
TEST_OUT_DIR=test/output

TEST_CASES=${wildcard ${TEST_CASES_DIR}/*.sh}
TEST_DIFF=${TEST_CASES:${TEST_CASES_DIR}/%.sh=${TEST_OUT_DIR}/%.diff}

.PHONY: clean distclean all

test: ${TEST_DIFF}

${TEST_OUT_DIR}/%.diff: ${TEST_CASES_DIR}/%.sh ${TEST_OUT_DIR}/%.dat
	@./$< | diff -B -w $(basename $@).dat - > $@; \
	if [ $$? -eq 0 ]; \
	then echo "$<" \
	     | awk '{ printf "%-40s \033[1m[OK]\033[m\n", $$1 }'; \
	else echo "$<" \
	     | awk '{ printf "%-40s \033[1m[ERROR]\033[m\n", $$1 }'; \
	fi

#
# Utilities
#

clean distclean:
	rm -f ${SHARED_OBJS} ${PERL_OBJS}
	rm -f ${TEST_DIFF}
