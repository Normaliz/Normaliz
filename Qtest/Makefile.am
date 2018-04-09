# TODO: Makefile.classic uses GNU make features.
# For portability it would be better to integrate Makefile.classic here and only use automake features.

TEST_TARGETS =
if ENFNORMALIZ
TEST_TARGETS += testnf
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
	test-nf/strict.ref

EXTRA_DIST += $(NMZTESTS)
