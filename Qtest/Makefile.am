# TODO: Makefile.classic uses GNU make features.
# For portability it would be better to integrate Makefile.classic here and only use automake features.

TEST_TARGETS = testqq testqqH
if ENFNORMALIZ
TEST_TARGETS += testnf testnfH
endif

check-local:
	rm -rf run_tests
	mkdir run_tests
	(cd $(srcdir) && tar cf - $(NMZTESTS)) | (cd run_tests && tar xf -)
	(cd run_tests && $(MAKE) -f $(abs_srcdir)/Makefile.classic QNORMALIZ=$(abs_top_builddir)/Qsource/Qnormaliz $(TEST_TARGETS))

clean-local:
	$(MAKE) -f $(srcdir)/Makefile.classic clean
	rm -rf run_tests

distclean-local:
	$(MAKE) -f $(srcdir)/Makefile.classic distclean
	rm -rf run_tests

EXTRA_DIST = Makefile.classic

NMZTESTS =					\
	test-nf/2cone.in			\
	test-nf/2cone.ref			\
	test-nf/2cone-sqrt2-cone.in		\
	test-nf/2cone-sqrt2-cone.ref		\
	test-nf/2cone-sqrt2.in			\
	test-nf/2cone-sqrt2.ref			\
	test-nf/600cell-dual.in			\
	test-nf/600cell-dual.ref		\
	test-nf/600cell.in			\
	test-nf/600cell.ref			\
	test-nf/buckyball-v.in			\
	test-nf/buckyball-v.ref			\
	test-nf/dodecahedron-v.in		\
	test-nf/dodecahedron-v.ref		\
	test-nf/formatted.in			\
	test-nf/formatted.ref			\
	test-nf/grandantiprism.in		\
	test-nf/grandantiprism.ref		\
	test-nf/great_rhombicuboctahedron-v.in	\
	test-nf/great_rhombicuboctahedron-v.ref	\
	test-nf/icosahedron-v.in		\
	test-nf/icosahedron-v.ref		\
	test-nf/ikosaeder32.in			\
	test-nf/ikosaeder32.ref			\
	test-nf/ikosaeder4.in			\
	test-nf/ikosaeder4.ref			\
	test-nf/ikosaeder.in			\
	test-nf/ikosaeder-perturbed-dual.in	\
	test-nf/ikosaeder-perturbed-dual.ref	\
	test-nf/ikosaeder-perturbed.in		\
	test-nf/ikosaeder-perturbed.ref		\
	test-nf/ikosaeder.ref			\
	test-nf/small_rhombicuboctahedron-v.in	\
	test-nf/small_rhombicuboctahedron-v.ref	\
	test-nf/sparse.in			\
	test-nf/sparse.ref			\
	test-nf/sqauare-3rd_root_2.in		\
	test-nf/sqauare-3rd_root_2.ref		\
	test-nf/sqauare-cone.in			\
	test-nf/sqauare-cone.ref		\
	test-nf/truncateddodecahedron.in	\
	test-nf/truncateddodecahedron.ref	\
	test-nf/dodecahedron-v.ref		\
	test-nf/formatted.in			\
	test-nf/formatted.ref			\
	test-nf/grandantiprism.in		\
	test-nf/grandantiprism.ref		\
	test-nf/great_rhombicuboctahedron-v.in	\
	test-nf/great_rhombicuboctahedron-v.ref	\
	test-nf/icosahedron-v.in		\
	test-nf/icosahedron-v.ref		\
	test-nf/ikosaeder32.in			\
	test-nf/ikosaeder32.ref			\
	test-nf/ikosaeder4.in			\
	test-nf/ikosaeder4.ref			\
	test-nf/ikosaeder.in			\
	test-nf/ikosaeder-perturbed-dual.in	\
	test-nf/ikosaeder-perturbed-dual.ref	\
	test-nf/ikosaeder-perturbed.in		\
	test-nf/ikosaeder-perturbed.ref		\
	test-nf/ikosaeder.ref			\
	test-nf/small_rhombicuboctahedron-v.in	\
	test-nf/small_rhombicuboctahedron-v.ref	\
	test-nf/QConstSymb.in			\
	test-nf/QConstSymb.ref			\
	test-nf/sparse.in			\
	test-nf/sparse.ref			\
	test-nf/sqauare-3rd_root_2.in		\
	test-nf/sqauare-3rd_root_2.ref		\
	test-nf/sqauare-cone.in			\
	test-nf/sqauare-cone.ref		\
	test-nf/truncateddodecahedron.in	\
	test-nf/truncateddodecahedron.ref	\
	test-nf/low-dim.in			\
	test-nf/low-dim.ref			\
	test-nf/nonpointed.in			\
	test-nf/nonpointed.ref			\
	test-nf/strict.in			\
	test-nf/strict.ref			\
	test-qq/2equations.in			\
	test-qq/2equations.ref			\
	test-qq/2x2.in				\
	test-qq/2x2.ref				\
	test-qq/3x3magic.in			\
	test-qq/3x3magic.ref			\
	test-qq/4x4.in				\
	test-qq/4x4.ref				\
	test-qq/5x5dual.in			\
	test-qq/5x5dual.ref			\
	test-qq/5x5.in				\
	test-qq/5x5.ref				\
	test-qq/A443.in				\
	test-qq/A443.ref			\
	test-qq/bo5.in				\
	test-qq/bo5.ref				\
	test-qq/Condorcet.in			\
	test-qq/Condorcet.ref			\
	test-qq/cut_poly_7.in			\
	test-qq/cut_poly_7.ref			\
	test-qq/gen_inhom_nonpointed.in		\
	test-qq/gen_inhom_nonpointed.ref	\
	test-qq/InhomIneq_7.in			\
	test-qq/InhomIneq_7.ref			\
	test-qq/InhomIneq_gen.in		\
	test-qq/InhomIneq_gen.ref		\
	test-qq/max_polytope_cand.in		\
	test-qq/max_polytope_cand.ref		\
	test-qq/medium.in			\
	test-qq/medium.ref			\
	test-qq/rational.in			\
	test-qq/rational.ref			\
	test-qq/small.in			\
	test-qq/small.ref			\
	test-qq/strict_sign.in			\
	test-qq/strict_sign.ref			\
	test-qq/subspace4.in			\
	test-qq/subspace4.ref
	
EXTRA_DIST += $(NMZTESTS)
