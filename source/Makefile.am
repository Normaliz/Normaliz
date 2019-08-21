AUTOMAKE_OPTIONS = subdir-objects std-options

EXTRA_DIST =

AM_CPPFLAGS = -I$(top_builddir)/source
AM_CXXFLAGS = -Wall -pedantic -Wno-unknown-pragmas -Wno-sign-compare
AM_LDFLAGS =

if ENABLE_OPENMP
AM_CXXFLAGS += $(OPENMP_CXXFLAGS)
endif
if NMZ_INTEGRATE
AM_CXXFLAGS += -DNMZ_COCOA -DCoCoA_THREADSAFE_HACK
endif

# Installed library, shared and/or static (using libtool).
lib_LTLIBRARIES = libnormaliz.la

# Installed headers
nobase_include_HEADERS = libnormaliz/cone.h libnormaliz/cone_property.h libnormaliz/convert.h libnormaliz/general.h libnormaliz/HilbertSeries.h libnormaliz/integer.h libnormaliz/input_type.h libnormaliz/matrix.h libnormaliz/my_omp.h libnormaliz/normaliz_exception.h libnormaliz/sublattice_representation.h libnormaliz/vector_operations.h libnormaliz/version.h libnormaliz/automorph.h  libnormaliz/libnormaliz.h  libnormaliz/map_operations.h libnormaliz/nmz_config.h
# Sources
libnormaliz_la_SOURCES = libnormaliz/enumeration.cpp libnormaliz/other_algorithms.cpp libnormaliz/linear_algebra.cpp libnormaliz/offload_handler.cpp libnormaliz/cone_and_control.cpp libnormaliz/primal.cpp libnormaliz/nmz_nauty.cpp libnormaliz/output.cpp
# Other headers (not installed)
noinst_HEADERS = libnormaliz/simplex.h libnormaliz/reduction.h libnormaliz/full_cone.h libnormaliz/cone_dual_mode.h libnormaliz/offload_handler.h libnormaliz/list_operations.h libnormaliz/bottom.h libnormaliz/project_and_lift.h libnormaliz/descent.h libnormaliz/nmz_nauty.h libnormaliz/nmz_integrate.h
# Sources included via libnormaliz-templated.cpp:
noinst_HEADERS += libnormaliz/input_type.cpp libnormaliz/integer.cpp libnormaliz/vector_operations.cpp libnormaliz/matrix.cpp libnormaliz/simplex.cpp libnormaliz/list_operations.cpp libnormaliz/sublattice_representation.cpp libnormaliz/reduction.cpp libnormaliz/full_cone.cpp libnormaliz/cone_dual_mode.cpp libnormaliz/cone.cpp libnormaliz/project_and_lift.cpp libnormaliz/HilbertSeries.cpp libnormaliz/cone_property.cpp libnormaliz/bottom.cpp libnormaliz/descent.cpp libnormaliz/output.cpp libnormaliz/general.cpp
noinst_HEADERS += libnormaliz/nmz_polynomial.cpp libnormaliz/nmz_integral.cpp libnormaliz/automorph.cpp libnormaliz/nmz_nauty.cpp

libnormaliz_la_LDFLAGS = -version-info 10:4:7
## Conjecture: x.y.z yields x-z.z.y
libnormaliz_la_LIBADD = $(GMP_LDFLAGS) -lgmpxx -lgmp

# Binary.
bin_PROGRAMS = normaliz

normaliz_SOURCES = normaliz.cpp
normaliz_LDADD = libnormaliz.la
noinst_HEADERS += options.h libnormaliz/output.h

# Sources included via Normaliz.cpp
noinst_HEADERS += input.cpp options.cpp

# Uninstalled binary maxsimplex
noinst_PROGRAMS = maxsimplex/maxsimplex
maxsimplex_maxsimplex_SOURCES = maxsimplex/maxsimplex.cpp
maxsimplex_maxsimplex_LDADD = libnormaliz.la

# Uninstalled binary dynamic
noinst_PROGRAMS += dynamic/dynamic
dynamic_dynamic_SOURCES = dynamic/dynamic.cpp
dynamic_dynamic_LDADD = libnormaliz.la

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

