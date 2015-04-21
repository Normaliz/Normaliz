## default normaliz bin
ifndef (NORMALIZ)
  NORMALIZ = ../source/normaliz
endif


#NORMPARA += -c
#NORMPARA += -x=1
DIFF = diff -b
#NICE = nice -n19


TEST_FILES_S = $(patsubst %.in,%.diff,$(wildcard test-s/*.in))
TEST_FILES_V = $(patsubst %.in,%.diff,$(wildcard test-v/*.in))
TEST_FILES_N = $(patsubst %.in,%.diff,$(wildcard test-n/*.in))
TEST_FILES_NN = $(patsubst %.in,%.diff,$(wildcard test-nn/*.in))
TEST_FILES_H = $(patsubst %.in,%.diff,$(wildcard test-h/*.in))
TEST_FILES_P = $(patsubst %.in,%.diff,$(wildcard test-p/*.in))
TEST_FILES_Q = $(patsubst %.in,%.diff,$(wildcard test-q/*.in))
TEST_FILES_1 = $(patsubst %.in,%.diff,$(wildcard test-1/*.in))
TEST_FILES_D = $(patsubst %.in,%.diff,$(wildcard test-d/*.in))
TEST_FILES_D1 = $(patsubst %.in,%.diff,$(wildcard test-d1/*.in))
TEST_FILES_R1 = $(patsubst %.in,%.diff,$(wildcard test-r1/*.in))
TEST_FILES_DEFAULT = $(patsubst %.in,%.diff,$(wildcard test-/*.in))


#fehlen noch: tests testr1
all: testv testn test1 testN testq testp testd testh testdefault testd1 testr1
tests: $(TEST_FILES_S)
testv: $(TEST_FILES_V)
testN: $(TEST_FILES_NN)
testn: $(TEST_FILES_N)
testh: $(TEST_FILES_H)
testp: $(TEST_FILES_P)
testq: $(TEST_FILES_Q)
test1: $(TEST_FILES_1)
testdefault: $(TEST_FILES_DEFAULT)
testd: $(TEST_FILES_D)
testd1: $(TEST_FILES_D1)
testr1: $(TEST_FILES_R1)

#compare rule
%.diff: %.out %.ref
	$(DIFF) $*.ref $*.out > $*.diff

#rules to generate the test files
test-s/%.out: $(NORMALIZ) test-s/%.in
	$(NICE) $(NORMALIZ) $(NORMPARA) -s test-s/$*
test-v/%.out: $(NORMALIZ) test-v/%.in
	$(NICE) $(NORMALIZ) $(NORMPARA) -v test-v/$*
test-1/%.out: $(NORMALIZ) test-1/%.in
	$(NICE) $(NORMALIZ) $(NORMPARA) -1 test-1/$*
test-n/%.out: $(NORMALIZ) test-n/%.in
	$(NICE) $(NORMALIZ) $(NORMPARA) -n test-n/$*
test-nn/%.out: $(NORMALIZ) test-nn/%.in
	$(NICE) $(NORMALIZ) $(NORMPARA) -N test-nn/$*
test-p/%.out: $(NORMALIZ) test-p/%.in
	$(NICE) $(NORMALIZ) $(NORMPARA) -p test-p/$*
test-q/%.out: $(NORMALIZ) test-q/%.in
	$(NICE) $(NORMALIZ) $(NORMPARA) -q test-q/$*
test-d/%.out: $(NORMALIZ) test-d/%.in
	$(NICE) $(NORMALIZ) $(NORMPARA) -d test-d/$*
test-h/%.out: $(NORMALIZ) test-h/%.in
	$(NICE) $(NORMALIZ) $(NORMPARA) -h test-h/$*
test-/%.out: $(NORMALIZ) test-/%.in
	$(NICE) $(NORMALIZ) $(NORMPARA) test-/$*
test-d1/%.out: $(NORMALIZ) test-d1/%.in
	$(NICE) $(NORMALIZ) $(NORMPARA) -d1 test-d1/$*
test-r1/%.out: $(NORMALIZ) test-r1/%.in
	$(NICE) $(NORMALIZ) $(NORMPARA) -r1 test-r1/$*



clean:
	-rm -f test-*/*.diff

distclean: clean
	-rm -f test-*/*.out

.PRECIOUS: test-s/%.out
.PRECIOUS: test-v/%.out
.PRECIOUS: test-1/%.out
.PRECIOUS: test-n/%.out
.PRECIOUS: test-nn/%.out
.PRECIOUS: test-p/%.out
.PRECIOUS: test-q/%.out
.PRECIOUS: test-d/%.out
.PRECIOUS: test-d1/%.out
.PRECIOUS: test-r1/%.out
.PRECIOUS: test-h/%.out
.PRECIOUS: test-/%.out
.PHONY: all testdefault tests testn test1 testN testp testh testd testd1 testH clean distclean
