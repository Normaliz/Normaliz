AUTOMAKE_OPTIONS = subdir-objects

EXTRA_DIST =

# Installed library, shared and/or static (using libtool).
lib_LTLIBRARIES = libnormaliz.la

# Installed headers
include_HEADERS = libnormaliz/cone.h libnormaliz/cone_property.h libnormaliz/convert.h libnormaliz/general.h libnormaliz/HilbertSeries.h libnormaliz/integer.h libnormaliz/libnormaliz.h libnormaliz/map_operations.h libnormaliz/matrix.h libnormaliz/my_omp.h libnormaliz/normaliz_exception.h libnormaliz/sublattice_representation.h libnormaliz/vector_operations.h libnormaliz/version.h
# Sources
libnormaliz_la_SOURCES = libnormaliz/cone_property.cpp libnormaliz/HilbertSeries.cpp libnormaliz/libnormaliz-templated.cpp libnormaliz/cone_helper.cpp libnormaliz/offload_handler.cpp libnormaliz/bottom.cpp 
# Other headers (not installed)
noinst_HEADERS = libnormaliz/simplex.h libnormaliz/reduction.h libnormaliz/full_cone.h libnormaliz/cone_dual_mode.h libnormaliz/offload_handler.h libnormaliz/list_operations.h libnormaliz/bottom.h libnormaliz/cone_helper.h
# Sources included via libnormaliz-templated.cpp:
noinst_HEADERS += libnormaliz/libnormaliz.cpp libnormaliz/integer.cpp libnormaliz/vector_operations.cpp libnormaliz/matrix.cpp libnormaliz/simplex.cpp libnormaliz/list_operations.cpp libnormaliz/sublattice_representation.cpp libnormaliz/reduction.cpp libnormaliz/full_cone.cpp libnormaliz/cone_dual_mode.cpp libnormaliz/cone.cpp

libnormaliz_la_LDFLAGS = -version-info 0:0:0
libnormaliz_la_LIBADD = $(GMP_LDFLAGS) -lgmpxx -lgmp

# Binary.
bin_PROGRAMS = normaliz

normaliz_SOURCES = Normaliz.cpp
normaliz_LDADD = libnormaliz.la
noinst_HEADERS += Normaliz.h Options.h output.h
# Sources included via Normaliz.cpp
noinst_HEADERS += Input.cpp Options.cpp output.cpp



## TODO: normaliz1 normaliz-pg

## Other things that go into the distribution:

# Cmake build system
EXTRA_DIST += CMakeLists.txt LibnormalizConfig.cmake.in UseLibnormaliz.cmake cmake/Modules cmake/detect_features.cmake cmake/Modules/FindGMP.cmake cmake/Modules/FindGMPXX.cmake cmake/Modules/FindMIC_GMP.cmake cmake/Modules/FindMIC_GMPXX.cmake cmake/Modules/FindSCIP.cmake libnormaliz/CMakeLists.txt libnormaliz/LibnormalizConfig.cmake

# Plain old build system
EXTRA_DIST += Makefile.classic Makefile.configuration libnormaliz/Makefile.classic
