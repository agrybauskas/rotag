#
# Unit tests
#

TEST_CASES_DIR=tests/cases
TEST_OUT_DIR=tests/outputs
TEST_CASES=${sort ${wildcard ${TEST_CASES_DIR}/*.sh}}
TEST_DIFF=${TEST_CASES:${TEST_CASES_DIR}/%.sh=${TEST_OUT_DIR}/%.diff}

.PHONY: clean distclean all

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
# Utilities
#

clean distclean:
	rm -f ${TEST_DIFF} ${VISUAL_TEST_JMOL}
