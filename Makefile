TEST_CASES_DIR=test/case
TEST_OUT_DIR=test/output

TEST_CASES=${wildcard ${TEST_CASES_DIR}/*.sh}
TEST_DIFF=${TEST_CASES:${TEST_CASES_DIR}/%.sh=${TEST_OUT_DIR}/%.diff}

.PHONY: clean distclean

test: ${TEST_DIFF} 

${TEST_OUT_DIR}/%.diff: ${TEST_CASES_DIR}/%.sh ${TEST_OUT_DIR}/%.dat
	@./$< | diff -B -w $(basename $@).dat -; \
	if [ $$? -eq 0 ]; \
	then echo $<"\t[OK]"; \
	else echo $<"\t[FALSE]"; \
	fi

clean distclean:
	rm -f ${TEST_DIFF}
