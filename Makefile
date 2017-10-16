#
# Compiling CPP and linking to Perl5 with SWIG.
#

CPP_DIR=lib/cpp
CPP_OBJ=${SWIG_FILE:%.i=%.o}
SWIG_FILE=${wildcard ${CPP_DIR}/*.i}
PM_FILE=${SWIG_FILE:%.i=%.pm}
WRAP_FILE=${SWIG_FILE:%.i=%_wrap.cxx}
WRAP_OBJ=${WRAP_FILE:%.cxx=%.o}
SHARED_OBJ=${CPP_OBJ:%.o=%.so}

build: ${SHARED_OBJ} ${PM_FILE}

%.pm %_wrap.cxx: %.i
	swig -c++ -perl $<

%.o: %.cpp
	g++ -c -fPIC $< -o $@

%_wrap.o: %_wrap.cxx
	g++ -c -fPIC $< -I$$(perl -e 'use Config; print $$Config{archlib};')/CORE \
	    -o $@

%.so: %_wrap.o %.o
	g++ -shared $^ -o $@

#
# Instalation of dependencies.
#

.PHONY: dependencies

dependencies:
	./dependencies/$$(lsb_release -is)-$$(lsb_release -rs)/install.sh

#
# Unit tests.
#

TEST_CASES_DIR=tests/cases
TEST_OUT_DIR=tests/outputs
TEST_CASES=${sort ${wildcard ${TEST_CASES_DIR}/*.sh}}
TEST_DIFF=${TEST_CASES:${TEST_CASES_DIR}/%.sh=${TEST_OUT_DIR}/%.diff}

test: ${TEST_DIFF}

${TEST_OUT_DIR}/%.diff: ${TEST_CASES_DIR}/%.sh ${TEST_OUT_DIR}/%.out
	@./$< 2>&1 | diff -B -w $(basename $@).out - > $@; \
	if [ $$? -eq 0 ]; \
	then echo "$<" \
	     | awk '{ printf "%-40s \033[1m[OK]\033[m\n",    $$1 }'; \
	else echo "$<" \
	     | awk '{ printf "%-40s \033[1m[ERROR]\033[m\n", $$1 }'; \
	       cat $@; \
	fi

#
# Utilities.
#

.PHONY: clean distclean

clean distclean:
	rm -f ${TEST_DIFF} ${VISUAL_TEST_JMOL} ${PM_FILE} ${SHARED_OBJ}
	rm -f ${CPP_OBJ} ${WRAP_FILE}
