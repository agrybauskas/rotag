#
# Unit tests
#

TEST_CASES_DIR=test/case
TEST_OUT_DIR=test/output
VISUAL_TEST_CASES_DIR=test/visual_case
VISUAL_TEST_OUT_DIR=test/output

TEST_CASES=${sort ${wildcard ${TEST_CASES_DIR}/*.sh}}
TEST_DIFF=${TEST_CASES:${TEST_CASES_DIR}/%.sh=${TEST_OUT_DIR}/%.diff}
VISUAL_TEST_CASES=${sort ${wildcard ${VISUAL_TEST_CASES_DIR}/*.sh}}
VISUAL_TEST_JMOL=${VISUAL_TEST_CASES:${VISUAL_TEST_CASES_DIR}/%.sh=${VISUAL_TEST_OUT_DIR}/%.jmol}

.PHONY: clean distclean all

test: ${TEST_DIFF}

visual_test: ${VISUAL_TEST_JMOL}

${TEST_OUT_DIR}/%.diff: ${TEST_CASES_DIR}/%.sh ${TEST_OUT_DIR}/%.dat
	@./$< 2>&1 | diff -B -w $(basename $@).dat - > $@; \
	if [ $$? -eq 0 ]; \
	then echo "$<" \
	     | awk '{ printf "%-40s \033[1m[OK]\033[m\n",    $$1 }'; \
	else echo "$<" \
	     | awk '{ printf "%-40s \033[1m[ERROR]\033[m\n", $$1 }'; \
	fi

${VISUAL_TEST_OUT_DIR}/%.jmol: ${VISUAL_TEST_CASES_DIR}/%.sh
	@./$< 2>&1 > $@; \
	if [ $$? -eq 0 ]; \
	then echo "$<" \
	     | awk '{ printf "%-40s \033[1m[DONE]\033[m\n",  $$1 }'; \
	else echo "$<" \
	     | awk '{ printf "%-40s \033[1m[ERROR]\033[m\n", $$1 }'; \
	fi

#
# Utilities
#

clean distclean:
	rm -f ${TEST_DIFF} ${VISUAL_TEST_JMOL}
