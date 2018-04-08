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
if ENFNORMALIZ
LIBS += -lnfemb -larb -lantic -lflint -lmpfr
endif
LIBS += -lgmpxx -lgmp

# Installed library, shared and/or static (using libtool).
lib_LTLIBRARIES = libQnormaliz.la

# Installed headers
nobase_include_HEADERS = libQnormaliz/Qcone.h libQnormaliz/Qcone_property.h libQnormaliz/Qconvert.h libQnormaliz/Qgeneral.h libQnormaliz/Qinteger.h libQnormaliz/libQnormaliz.h libQnormaliz/Qmap_operations.h libQnormaliz/Qmatrix.h libQnormaliz/Qmy_omp.h libQnormaliz/Qnormaliz_exception.h libQnormaliz/Qsublattice_representation.h libQnormaliz/Qvector_operations.h libQnormaliz/Qversion.h
# Sources
libQnormaliz_la_SOURCES = libQnormaliz/Qcone_property.cpp libQnormaliz/libQnormaliz-templated.cpp libQnormaliz/Qcone_helper.cpp 
# Other headers (not installed)
noinst_HEADERS = libQnormaliz/Qfull_cone.h libQnormaliz/Qlist_operations.h libQnormaliz/Qcone_helper.h
# Sources included via libQnormaliz-templated.cpp:
noinst_HEADERS += libQnormaliz/Qinteger.cpp libQnormaliz/Qvector_operations.cpp libQnormaliz/Qmatrix.cpp libQnormaliz/Qlist_operations.cpp libQnormaliz/Qsublattice_representation.cpp libQnormaliz/Qfull_cone.cpp libQnormaliz/Qcone.cpp libQnormaliz/Qlibnormaliz.cpp

libQnormaliz_la_LDFLAGS = $(AM_LDFLAGS) -version-info 6:0:3
## yields version number 3.2.0 for libQnormaliz.so on Linux
## Conjecture: x.y.z yields x-z.z.y
libQnormaliz_la_LIBADD = ../source/libnormaliz.la $(GMP_LDFLAGS) -lgmpxx -lgmp

# Binary.
bin_PROGRAMS = Qnormaliz

Qnormaliz_SOURCES = Qnormaliz.cpp
Qnormaliz_LDADD = libQnormaliz.la
noinst_HEADERS += Qnormaliz.h Qoptions.h Qoutput.h
# Sources included via Qnormaliz.cpp
noinst_HEADERS += Qinput.cpp Qoptions.cpp Qoutput.cpp

## Other things that go into the distribution:

# Cmake build system
EXTRA_DIST += CMakeLists.txt LibQnormalizConfig.cmake.in UseLibQnormaliz.cmake cmake/Modules cmake/detect_features.cmake cmake/Modules/FindGMP.cmake cmake/Modules/FindGMPXX.cmake cmake/Modules/FindMIC_GMP.cmake cmake/Modules/FindMIC_GMPXX.cmake cmake/Modules/FindSCIP.cmake libQnormaliz/CMakeLists.txt

# Plain old build system
EXTRA_DIST += Makefile.classic Makefile.configuration libQnormaliz/Makefile.classic

