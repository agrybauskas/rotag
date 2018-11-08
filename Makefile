#
# Grammar compilation.
#

YAPP_DIR=lib/Grammar
YAPP_FILES=${wildcard ${YAPP_DIR}/*.yp}
GRAMMAR_MODULES=${YAPP_FILES:%.yp=%.pm}

all: ${GRAMMAR_MODULES} ${PERL_MODULE}

%.pm: %.yp
	yapp -o $@ $<

#
# Generate Perl modules.
#

LIB_DIR=lib
TOOLS_DIR=tools
PERL_TEMPLATE=${LIB_DIR}/Constants.pmin
PERL_MODULE=${LIB_DIR}/Constants.pm

all: ${PERL_MODULE}

${PERL_MODULE}: ${PERL_TEMPLATE}
	sed 's/@PI@/'$$(${TOOLS_DIR}/calculate-pi)'/g' $^ \
	    | sed 's/@EPSILON@/'$$(${TOOLS_DIR}/calculate-epsilon)'/g' \
	    > $@

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

test: ${GRAMMAR_MODULES} ${PERL_MODULE} | ${TEST_DIFF}

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

clean:
	rm -f ${TEST_DIFF}
	rm -f ${COVERAGE_CASES}
	rm -f ${COVERAGE_OUTS}
	rm -fr ${COVERAGE_CASES_DIR}/cover_db

cleanAll distclean: clean
	rm -f ${GRAMMAR_MODULES}
	rm -f ${PERL_MODULE}
