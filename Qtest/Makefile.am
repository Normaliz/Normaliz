# TODO: Makefile.classic uses GNU make features.
# For portability it would be better to integrate Makefile.classic here and only use automake features.

check-local:
	rm -rf run_tests
	mkdir run_tests
	(cd $(srcdir) && tar cf - $(NMZTESTS)) | (cd run_tests && tar xf -)
	(cd run_tests && $(MAKE) -f $(abs_srcdir)/Makefile.classic QNORMALIZ=$(abs_top_builddir)/Qsource/Qnormaliz)

clean-local:
	$(MAKE) -f $(srcdir)/Makefile.classic clean
	rm -rf run_tests

distclean-local:
	$(MAKE) -f $(srcdir)/Makefile.classic distclean
	rm -rf run_tests

EXTRA_DIST = Makefile.classic

NMZTESTS =					\
	test-s/2equations.in			\
	test-s/2equations.ref			\
	test-s/2x2.in				\
	test-s/2x2.ref				\
	test-s/3x3magic.in			\
	test-s/3x3magic.ref			\
	test-s/4x4.in				\
	test-s/4x4.ref				\
	test-s/5x5.in				\
	test-s/5x5.ref				\
	test-s/5x5dual.in			\
	test-s/5x5dual.ref			\
	test-s/A443.in				\
	test-s/A443.ref				\
	test-s/bo5.in				\
	test-s/bo5.ref				\
	test-s/Condorcet.in			\
	test-s/Condorcet.ref			\
	test-s/cut_poly_7.in			\
	test-s/cut_poly_7.ref			\
	test-s/gen_inhom_nonpointed.in		\
	test-s/gen_inhom_nonpointed.ref		\
	test-s/InhomIneq_7.in			\
	test-s/InhomIneq_7.ref			\
	test-s/InhomIneq_gen.in			\
	test-s/InhomIneq_gen.ref		\
	test-s/max_polytope_cand.in		\
	test-s/max_polytope_cand.ref		\
	test-s/medium.in			\
	test-s/medium.ref			\
	test-s/rational.in			\
	test-s/rational.ref			\
	test-s/small.in				\
	test-s/small.ref			\
	test-s/strict_sign.in			\
	test-s/strict_sign.ref			\
	test-s/subspace4.in			\
	test-s/subspace4.ref

EXTRA_DIST += $(NMZTESTS)
