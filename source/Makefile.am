AUTOMAKE_OPTIONS = subdir-objects std-options

EXTRA_DIST =

AM_CPPFLAGS = -I$(top_builddir)/source
AM_CXXFLAGS = 
AM_LDFLAGS =

if ENABLE_OPENMP
AM_CXXFLAGS += $(OPENMP_CXXFLAGS)
endif
if NMZ_SCIP
AM_CXXFLAGS += $(SCIP_CXXFLAGS)
AM_LDFLAGS += $(SCIP_LDFLAGS)
LIBS += $(SCIP_LIBS)
endif
if NMZ_INTEGRATE
AM_CXXFLAGS += -DNMZ_COCOA -DCoCoA_THREADSAFE_HACK
LIBS += $(COCOA_LIBS)
endif
if NMZ_FLINT
LIBS += $(FLINT_LIBS)

endif
LIBS += -lgmpxx -lgmp -lpthread

# Installed library, shared and/or static (using libtool).
lib_LTLIBRARIES = libnormaliz.la

# Installed headers
nobase_include_HEADERS = libnormaliz/cone.h libnormaliz/cone_property.h libnormaliz/convert.h libnormaliz/general.h libnormaliz/HilbertSeries.h libnormaliz/integer.h libnormaliz/libnormaliz.h libnormaliz/map_operations.h libnormaliz/matrix.h libnormaliz/my_omp.h libnormaliz/normaliz_exception.h libnormaliz/sublattice_representation.h libnormaliz/vector_operations.h libnormaliz/version.h libnormaliz/nmz_integrate.h
# Sources
libnormaliz_la_SOURCES = libnormaliz/enumeration.cpp libnormaliz/other_algorithms.cpp libnormaliz/linear_algebra.cpp libnormaliz/offload_handler.cpp libnormaliz/cone_and_control.cpp libnormaliz/primal.cpp
# Other headers (not installed)
noinst_HEADERS = libnormaliz/simplex.h libnormaliz/reduction.h libnormaliz/full_cone.h libnormaliz/cone_dual_mode.h libnormaliz/offload_handler.h libnormaliz/list_operations.h libnormaliz/bottom.h libnormaliz/project_and_lift.h libnormaliz/descent.h
# Sources included via libnormaliz-templated.cpp:
noinst_HEADERS += libnormaliz/libnormaliz.cpp libnormaliz/integer.cpp libnormaliz/vector_operations.cpp libnormaliz/matrix.cpp libnormaliz/simplex.cpp libnormaliz/list_operations.cpp libnormaliz/sublattice_representation.cpp libnormaliz/reduction.cpp libnormaliz/full_cone.cpp libnormaliz/cone_dual_mode.cpp libnormaliz/cone.cpp libnormaliz/project_and_lift.cpp libnormaliz/HilbertSeries.cpp libnormaliz/nmz_integrate.cpp libnormaliz/cone_property.cpp libnormaliz/bottom.cpp libnormaliz/descent.cpp
# Sources included via nmz_integrate.cpp:
noinst_HEADERS += libnormaliz/nmz_polynomial.cpp libnormaliz/nmz_integral.cpp

libnormaliz_la_LDFLAGS = -version-info 9:1:6
## yields version number 3.6.1 for libnormaliz.so on Linux
## Conjecture: x.y.z yields x-z.z.y
libnormaliz_la_LIBADD = $(GMP_LDFLAGS) -lgmpxx -lgmp

# Binary.
bin_PROGRAMS = normaliz

normaliz_SOURCES = normaliz.cpp
normaliz_LDADD = libnormaliz.la
noinst_HEADERS += normaliz.h options.h output.h
# Sources included via Normaliz.cpp
noinst_HEADERS += input.cpp options.cpp output.cpp

# Uninstalled binary maxsimplex
noinst_PROGRAMS = maxsimplex/maxsimplex
maxsimplex_maxsimplex_SOURCES = maxsimplex/maxsimplex.cpp
maxsimplex_maxsimplex_LDADD = libnormaliz.la

## Other things that go into the distribution:

# Cmake build system
EXTRA_DIST += CMakeLists.txt LibnormalizConfig.cmake.in UseLibnormaliz.cmake cmake/Modules cmake/detect_features.cmake cmake/Modules/FindGMP.cmake cmake/Modules/FindGMPXX.cmake cmake/Modules/FindMIC_GMP.cmake cmake/Modules/FindMIC_GMPXX.cmake cmake/Modules/FindSCIP.cmake libnormaliz/CMakeLists.txt libnormaliz/LibnormalizConfig.cmake maxsimplex/CMakeLists.txt

# Plain old build system
EXTRA_DIST += Makefile.classic Makefile.configuration libnormaliz/Makefile.classic maxsimplex/Makefile.classic

# The classic makefile creates this file, which can confuse automake builds
CLEANFILES = normaliz.o

# Documentation

EXTRA_DIST += libnormaliz/README

EXTRA_DIST += outerpar/outerpar.cpp
EXTRA_DIST += outerpar/Makefile.classic
EXTRA_DIST += outerpar/small_gens.mat

