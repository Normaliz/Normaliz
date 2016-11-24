AUTOMAKE_OPTIONS = subdir-objects

EXTRA_DIST =

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
LIBS += -lgmpxx -lgmp

# Installed library, shared and/or static (using libtool).
lib_LTLIBRARIES = libnormaliz.la

# Installed headers
nobase_include_HEADERS = libnormaliz/cone.h libnormaliz/cone_property.h libnormaliz/convert.h libnormaliz/general.h libnormaliz/HilbertSeries.h libnormaliz/integer.h libnormaliz/libnormaliz.h libnormaliz/map_operations.h libnormaliz/matrix.h libnormaliz/my_omp.h libnormaliz/normaliz_exception.h libnormaliz/sublattice_representation.h libnormaliz/vector_operations.h libnormaliz/version.h
# Sources
libnormaliz_la_SOURCES = libnormaliz/cone_property.cpp libnormaliz/HilbertSeries.cpp libnormaliz/libnormaliz-templated.cpp libnormaliz/cone_helper.cpp libnormaliz/offload_handler.cpp libnormaliz/bottom.cpp 
# Other headers (not installed)
noinst_HEADERS = libnormaliz/simplex.h libnormaliz/reduction.h libnormaliz/full_cone.h libnormaliz/cone_dual_mode.h libnormaliz/offload_handler.h libnormaliz/list_operations.h libnormaliz/bottom.h libnormaliz/cone_helper.h
# Sources included via libnormaliz-templated.cpp:
noinst_HEADERS += libnormaliz/libnormaliz.cpp libnormaliz/integer.cpp libnormaliz/vector_operations.cpp libnormaliz/matrix.cpp libnormaliz/simplex.cpp libnormaliz/list_operations.cpp libnormaliz/sublattice_representation.cpp libnormaliz/reduction.cpp libnormaliz/full_cone.cpp libnormaliz/cone_dual_mode.cpp libnormaliz/cone.cpp

libnormaliz_la_LDFLAGS = -version-info 4:4:1 
## yields version number 3.1.2 for libnormaliz.so on Linux
libnormaliz_la_LIBADD = $(GMP_LDFLAGS) -lgmpxx -lgmp

# Binary.
bin_PROGRAMS = normaliz

normaliz_SOURCES = normaliz.cpp
normaliz_LDADD = libnormaliz.la
noinst_HEADERS += normaliz.h options.h output.h
# Sources included via Normaliz.cpp
noinst_HEADERS += input.cpp options.cpp output.cpp

# nmzIntegrate.
if NMZ_INTEGRATE
bin_PROGRAMS += nmzIntegrate
nmzIntegrate_SOURCES = genEhrhart/nmzIntegrate.C
nmzIntegrate_CPPFLAGS = -IgenEhrhart $(AM_CPPFLAGS)
# included in nmzIntegrate.C
noinst_HEADERS += genEhrhart/nmzIntInput.C genEhrhart/nmzIntPoly.C genEhrhart/cyclRatFunct.C genEhrhart/genEhrhart.C genEhrhart/nmzIntegral.C
nmzIntegrate_LDADD = $(COCOA_LIBS) libnormaliz.la
endif

# Uninstalled binary maxsimplex
noinst_PROGRAMS = maxsimplex/maxsimplex
maxsimplex_maxsimplex_SOURCES = maxsimplex/maxsimplex.cpp
maxsimplex_maxsimplex_LDADD = libnormaliz.la

## Other things that go into the distribution:

# Cmake build system
EXTRA_DIST += CMakeLists.txt LibnormalizConfig.cmake.in UseLibnormaliz.cmake cmake/Modules cmake/detect_features.cmake cmake/Modules/FindGMP.cmake cmake/Modules/FindGMPXX.cmake cmake/Modules/FindMIC_GMP.cmake cmake/Modules/FindMIC_GMPXX.cmake cmake/Modules/FindSCIP.cmake libnormaliz/CMakeLists.txt libnormaliz/LibnormalizConfig.cmake

# Plain old build system
EXTRA_DIST += Makefile.classic Makefile.configuration libnormaliz/Makefile.classic

# The classic makefile creates this file, which can confuse automake builds
CLEANFILES = normaliz.o

# Documentation

# TODO: shouldn't these two be moved to toplevel?
EXTRA_DIST += COPYING INSTALL

EXTRA_DIST += libnormaliz/README

