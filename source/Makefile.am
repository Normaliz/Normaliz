EXTRA_DIST =

AM_CPPFLAGS = -I$(top_srcdir)/source
AM_CXXFLAGS = -Wall -pedantic -Wno-unknown-pragmas -Wno-sign-compare
AM_CXXFLAGS += $(OPENMP_CXXFLAGS)
AM_LDFLAGS = -Wl,-s


LIBS += $(COCOA_LIBS)
LIBS += $(E_ANTIC_LIBS)
LIBS += $(FLINT_LIBS)
LIBS += $(NAUTY_LIBS)
LIBS += $(GMP_LIBS)


# Installed library, shared and/or static (using libtool).
lib_LTLIBRARIES = libnormaliz.la

# Installed headers
nobase_include_HEADERS = \
  libnormaliz/automorph.h \
  libnormaliz/cone_property.h \
  libnormaliz/cone.h \
  libnormaliz/descent.h \
  libnormaliz/dynamic_bitset.h \
  libnormaliz/fusion.h \
  libnormaliz/general.h \
  libnormaliz/HilbertSeries.h \
  libnormaliz/input.h \
  libnormaliz/input_type.h \
  libnormaliz/integer.h \
  libnormaliz/libnormaliz.h \
  libnormaliz/list_and_map_operations.h \
  libnormaliz/matrix.h \
  libnormaliz/my_omp.h \
  libnormaliz/nmz_config.h \
  libnormaliz/nmz_integrate.h \
  libnormaliz/nmz_nauty.h \
  libnormaliz/nmz_polynomial.h \
  libnormaliz/normaliz_exception.h \
  libnormaliz/options.h \
  libnormaliz/output.h \
  libnormaliz/sublattice_representation.h \
  libnormaliz/vector_operations.h \
  libnormaliz/version.h

# Sources
libnormaliz_la_SOURCES = \
  libnormaliz/automorph.cpp \
  libnormaliz/binomial.cpp \
  libnormaliz/binomial_containers.cpp \
  libnormaliz/chunk.cpp \
  libnormaliz/collection.cpp \
  libnormaliz/collect_lat.cpp \
  libnormaliz/cone_dual_mode.cpp \
  libnormaliz/cone_property.cpp \
  libnormaliz/cone.cpp \
  libnormaliz/descent.cpp \
  libnormaliz/face_lattice.cpp \
  libnormaliz/full_cone.cpp \
  libnormaliz/fusion.cpp \
  libnormaliz/general.cpp \
  libnormaliz/HilbertSeries.cpp \
  libnormaliz/induction.cpp \
  libnormaliz/input.cpp \
  libnormaliz/lattice_ideal.cpp \
  libnormaliz/matrix.cpp \
  libnormaliz/nmz_integral.cpp \
  libnormaliz/nmz_nauty.cpp \
  libnormaliz/nmz_polynomial.cpp \
  libnormaliz/nmz_hash.cpp \
  libnormaliz/offload_handler.cpp \
  libnormaliz/options.cpp \
  libnormaliz/output.cpp \
  libnormaliz/project_and_lift.cpp \
  libnormaliz/reduction.cpp \
  libnormaliz/signed_dec.cpp \
  libnormaliz/simplex.cpp \
  libnormaliz/sublattice_representation.cpp


# Other headers (not installed)
noinst_HEADERS = \
  libnormaliz/binomial.h \
  libnormaliz/binomial_containers.h \
  libnormaliz/chunk.h \
  libnormaliz/collect_lat.h \
  libnormaliz/cone_dual_mode.h \
  libnormaliz/face_lattice.h \
  libnormaliz/full_cone.h \
  libnormaliz/induction.h \
  libnormaliz/lattice_ideal.h \
  libnormaliz/nmz_hash.h \
  libnormaliz/offload_handler.h \
  libnormaliz/options.h \
  libnormaliz/project_and_lift.h \
  libnormaliz/reduction.h \
  libnormaliz/signed_dec.h \
  libnormaliz/simplex.h \
  libnormaliz/collection.h

libnormaliz_la_LDFLAGS = -no-undefined -version-info 13:3:10
## Conjecture: x.y.z yields x-z.z.y
#libnormaliz_la_LIBADD = $(GMP_LIBS)

# Binary.
bin_PROGRAMS = normaliz

normaliz_SOURCES = normaliz.cpp
normaliz_LDADD = libnormaliz.la

# Uninstalled binary maxsimplex
noinst_PROGRAMS = maxsimplex/maxsimplex
maxsimplex_maxsimplex_SOURCES = maxsimplex/maxsimplex.cpp
maxsimplex_maxsimplex_LDADD = libnormaliz.la

# Uninstalled binary dynamic
noinst_PROGRAMS += dynamic/dynamic
dynamic_dynamic_SOURCES = dynamic/dynamic.cpp
dynamic_dynamic_LDADD = libnormaliz.la

EXTRA_DIST += dynamic/dynamic.cpp
EXTRA_DIST += dynamic/second.mat
EXTRA_DIST += dynamic/third.mat

## Other things that go into the distribution:


# Plain old build system
EXTRA_DIST += Makefile.classic Makefile.configuration libnormaliz/Makefile.classic

# The classic makefile creates this file, which can confuse automake builds
CLEANFILES = normaliz.o

# Documentation

EXTRA_DIST += libnormaliz/README

EXTRA_DIST += outerpar/outerpar.cpp
EXTRA_DIST += outerpar/small_gens.mat

