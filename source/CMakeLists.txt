# Pre-requisites about cmake itself
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0005 NEW) # new in 2.6 
  cmake_policy(SET CMP0003 NEW) # new in 2.6 
  if(POLICY CMP0042)            # new in 2.8.2
    cmake_policy(SET CMP0042 NEW)  # MACOSX_RPATH enabled by default
  endif(POLICY CMP0042)
endif(COMMAND cmake_policy)
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required

# The project name and the type of project
PROJECT(libnormaliz)
SET(VERSION_MAJOR 3)
SET(VERSION_MINOR 5)
SET(VERSION_PATCH 1)
SET(VERSION_STRING ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})
SET(VERSION ${VERSION_STRING})

# automatically set lib suffix, disabled for now
#IF(UNIX AND NOT APPLE AND NOT CMAKE_CROSSCOMPILING AND NOT EXISTS "/etc/debian_version")
#  IF(CMAKE_SIZEOF_VOID_P EQUAL 8 AND NOT LIB_SUFFIX)
#    SET(LIB_SUFFIX 64)
#  ENDIF()
#ENDIF()

# PkgConfig
SET(prefix      ${CMAKE_INSTALL_PREFIX})
SET(exec_prefix ${CMAKE_INSTALL_PREFIX}/bin)
SET(includedir  ${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME})

IF(${UNIX})
   SET(libdir      ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX})
ENDIF(${UNIX})
IF(${WIN32})
   SET(libdir      ${CMAKE_INSTALL_PREFIX}/bin)
ENDIF(${WIN32})

#set default build type
#IF (NOT CMAKE_BUILD_TYPE)  #or IF (NOT DEFINED CMAKE_BUILD_TYPE)
#  SET(CMAKE_BUILD_TYPE Release)
#ENDIF (NOT CMAKE_BUILD_TYPE)

##### Always full RPATH  #####

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

##############################

SET(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib") #needed on mac


# compiler flags
IF (MSVC)
  SET (CXX_WARNING_FLAGS "/Wall")
  SET (CXX_OPTIMIZATION_FLAGS "/O2")
ELSE ()
  SET (CXX_WARNING_FLAGS "-Wall -pedantic")
  SET (CXX_OPTIMIZATION_FLAGS "-O3")
ENDIF ()

IF (NOT MSVC)
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
  endif()
  IF (${APPLE})
    IF (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND (COMPILER_SUPPORTS_CXX11 OR COMPILER_SUPPORTS_CXX0X) )
      # using clang with new c++11 library on Apple
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    ENDIF ()
  ENDIF ()
ENDIF()


SET (NMZ_OPENMP "YES" CACHE BOOL "use openmp")
IF (NMZ_OPENMP)
  FIND_PACKAGE(OpenMP)
  IF (NOT OPENMP_FOUND)
    MESSAGE(STATUS "Not using OpenMP, compiler does not support it")
    SET (NMZ_OPENMP "NO" CACHE BOOL "use openmp" FORCE)
  ENDIF(NOT OPENMP_FOUND)
ENDIF (NMZ_OPENMP)
IF (NMZ_OPENMP)
  include("${CMAKE_SOURCE_DIR}/cmake/detect_features.cmake")
  NMZ_CHECK_FOR_EXCEPTION_PTR()
  IF (NOT HAVE_EXCEPTION_PTR)
    MESSAGE(STATUS "Not using OpenMP, compiler does not support std::exception_ptr and std::rethrow")
    SET (NMZ_OPENMP "NO" CACHE BOOL "use openmp" FORCE)
  ENDIF (NOT HAVE_EXCEPTION_PTR)
ENDIF (NMZ_OPENMP)
IF (NMZ_OPENMP)
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ELSE (NMZ_OPENMP)
  SET (CXX_WARNING_FLAGS "${CXX_WARNING_FLAGS} -Wno-unknown-pragmas")
ENDIF (NMZ_OPENMP)

# add compiler flags
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_WARNING_FLAGS} ${CXX_OPTIMIZATION_FLAGS}")

INCLUDE_DIRECTORIES(.)

# handle static compilation
SET (NMZ_BUILD_STATIC "NO" CACHE BOOL "build normaliz statically")
IF (NMZ_BUILD_STATIC)
  IF (${APPLE})
    IF (CMAKE_COMPILER_IS_GNUCXX)
        SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
    ENDIF (CMAKE_COMPILER_IS_GNUCXX)
    # make it compatible with older Mac versions
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmacosx-version-min=10.7")
  ELSEIF (${UNIX})
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  ELSE ()
    MESSAGE (FATAL_ERROR "Static linking not (yet) supported on this system")
  ENDIF ()
ENDIF (NMZ_BUILD_STATIC)


# handle usage of SCIP
SET (NMZ_SCIP "YES" CACHE BOOL "Use SCIP")
FIND_PACKAGE(SCIP) # should always be defined
IF(SCIP_FOUND)
  IF(NMZ_SCIP)
    SET (NMZ_SCIP_FLAGS "-DNMZ_SCIP")
    INCLUDE_DIRECTORIES(${SCIP_INCLUDE_DIR})
    MESSAGE(STATUS "Adding to compiler flags NMZ_SCIP_FLAGS=${NMZ_SCIP_FLAGS}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NMZ_SCIP_FLAGS}")
  ENDIF(NMZ_SCIP)
ELSE(SCIP_FOUND)
    SET (SCIP_LIBRARY "")
ENDIF(SCIP_FOUND)


# handle use of CoCoALib
SET (NMZ_COCOA "YES" CACHE BOOL "Use CoCoALib")
FIND_PACKAGE(COCOA)
IF(COCOA_FOUND)
  IF(NMZ_COCOA)
    SET (NMZ_COCOA_FLAGS "-DNMZ_COCOA -DCoCoA_THREADSAFE_HACK")
    INCLUDE_DIRECTORIES(${COCOA_INCLUDE_DIR})
    MESSAGE(STATUS "Adding to compiler flags NMZ_COCOA_FLAGS=${NMZ_COCOA_FLAGS}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NMZ_COCOA_FLAGS}")
  ENDIF(NMZ_COCOA)
ELSE(COCOA_FOUND)
    SET (COCOA_LIBRARY "")
ENDIF (COCOA_FOUND)

# handle use of Flint
SET (NMZ_FLINT "YES" CACHE BOOL "Use Flint")
FIND_PACKAGE(MPFR)
FIND_PACKAGE(FLINT)
IF(FLINT_FOUND AND MPFR_FOUND)
  IF(NMZ_FLINT)
    SET (NMZ_FLINT_FLAGS "-DNMZ_FLINT ")
    INCLUDE_DIRECTORIES(${FLINT_INCLUDE_DIR})
    INCLUDE_DIRECTORIES($${MPFR_INCLUDE_DIR})
    MESSAGE(STATUS "Adding to compiler flags NMZ_FLINT_FLAGS=${NMZ_FLINT_FLAGS}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NMZ_FLINT_FLAGS}")
  ENDIF(NMZ_FLINT)
ELSE(FLINT_FOUND AND MPFR_FOUND)
    SET (FLINT_LIBRARY "")
    SET(MPFR_LIBRARY " ")
ENDIF (FLINT_FOUND AND MPFR_FOUND)

# find gmp
FIND_PACKAGE(GMP REQUIRED)
INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})
IF (GMP_FOUND)
  MESSAGE(STATUS "GMP found with include dir: ${GMP_INCLUDE_DIR}")
  MESSAGE(STATUS "GMP found with lib dir: ${GMP_LIBRARIES_DIR}")
  MESSAGE(STATUS "GMP found with lib : ${GMP_LIBRARIES}")
  MESSAGE(STATUS "GMP found with static lib : ${GMP_STATIC_LIBRARIES}")
ENDIF()

FIND_PACKAGE(GMPXX REQUIRED)
INCLUDE_DIRECTORIES(${GMPXX_INCLUDE_DIR})
IF (GMPXX_FOUND)
  MESSAGE(STATUS "GMPXX found with include dir: ${GMPXX_INCLUDE_DIR}")
  MESSAGE(STATUS "GMPXX found with lib : ${GMPXX_LIBRARIES}")
  MESSAGE(STATUS "GMPXX found with static lib : ${GMPXX_STATIC_LIBRARIES}")
ENDIF()


# handle offload to mic (intel xeon phi)
SET (NMZ_MIC_OFFLOAD "NO" CACHE BOOL "Use offload to mic (intel xeon phi)")
IF (NMZ_MIC_OFFLOAD)
  IF ( NOT UNIX OR APPLE )
    MESSAGE (FATAL_ERROR "Offload to mic only supported on Linux")
  ENDIF()
  IF ( NMZ_SCIP AND SCIP_FOUND ) 
    MESSAGE (FATAL_ERROR "Simultaneous use of offload and scip not supported")
  ENDIF()

  # find gmp for mic
  FIND_PACKAGE(MIC_GMP REQUIRED)
  INCLUDE_DIRECTORIES(${MIC_GMP_INCLUDE_DIR})  #probably not the correct place
  IF (MIC_GMP_FOUND)
    MESSAGE(STATUS "MIC_GMP found with include dir: ${MIC_GMP_INCLUDE_DIR}")
    MESSAGE(STATUS "MIC_GMP found with lib dir: ${MIC_GMP_LIBRARIES_DIR}")
    MESSAGE(STATUS "MIC_GMP found with lib : ${MIC_GMP_LIBRARIES}")
    MESSAGE(STATUS "MIC_GMP found with static lib : ${MIC_GMP_STATIC_LIBRARIES}")
  ENDIF()
  FIND_PACKAGE(MIC_GMPXX REQUIRED)
  INCLUDE_DIRECTORIES(${MIC_GMPXX_INCLUDE_DIR}) #probably not the correct place
  IF (MIC_GMPXX_FOUND)
    MESSAGE(STATUS "MIC_GMPXX found with include dir: ${MIC_GMPXX_INCLUDE_DIR}")
    MESSAGE(STATUS "MIC_GMPXX found with lib : ${MIC_GMPXX_LIBRARIES}")
    MESSAGE(STATUS "MIC_GMPXX found with static lib : ${MIC_GMPXX_STATIC_LIBRARIES}")
  ENDIF()
  SET (NMZ_MIC_FLAGS "-offload-option,mic,compiler,\"-L${MIC_GMP_LIBRARIES_DIR} -lgmpxx -lgmp\" -DNMZ_MIC_OFFLOAD")
  SET (NMZ_MIC_FLAGS "${NMZ_MIC_FLAGS} -wd3218")
  MESSAGE(STATUS "Adding to compiler flags NMZ_MIC_FLAGS=${NMZ_MIC_FLAGS}")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NMZ_MIC_FLAGS}")
ENDIF (NMZ_MIC_OFFLOAD)


ADD_SUBDIRECTORY(libnormaliz)

# add the normaliz executable
ADD_EXECUTABLE (normaliz normaliz.cpp)

IF (NMZ_BUILD_STATIC)
  TARGET_LINK_LIBRARIES (normaliz normaliz-static)
  # TARGET_LINK_LIBRARIES (normaliz ${GMPXX_STATIC_LIBRARIES} ${GMP_STATIC_LIBRARIES})
  #Static start and end
  SET_TARGET_PROPERTIES(normaliz PROPERTIES LINK_SEARCH_START_STATIC 1)
  SET_TARGET_PROPERTIES(normaliz PROPERTIES LINK_SEARCH_END_STATIC 1)
  # Don't try to modify the RPATH
  SET_TARGET_PROPERTIES(normaliz PROPERTIES INSTALL_RPATH "")
ELSE (NMZ_BUILD_STATIC)
  TARGET_LINK_LIBRARIES (normaliz normaliz-shared)
  # TARGET_LINK_LIBRARIES (normaliz ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
ENDIF (NMZ_BUILD_STATIC)

INSTALL(TARGETS normaliz DESTINATION bin)

SET (NMZ_MAX "YES" CACHE BOOL "Compile with maxsimplex")
IF (NMZ_MAX)
    ADD_SUBDIRECTORY(maxsimplex)
    #TARGET_LINK_LIBRARIES (normaliz maxsimplex) 
ENDIF (NMZ_MAX)

#custom make targets
ADD_CUSTOM_TARGET (uninstall @echo uninstall package)
IF (UNIX)
  ADD_CUSTOM_COMMAND(
    COMMENT "uninstall package"
    COMMAND xargs ARGS rm < install_manifest.txt
    TARGET  uninstall
  )
ELSE(UNIX)
  ADD_CUSTOM_COMMAND(
    COMMENT "uninstall only implemented in unix"
    TARGET  uninstall
  )
ENDIF(UNIX)

ADD_CUSTOM_TARGET (dep @echo install package dependencies)
IF (UNIX)
  ADD_CUSTOM_COMMAND(
    COMMENT "make package dependencies"
    COMMAND ./make_dep.sh ARGS -d ../dependencies
    TARGET  dep
  )
ELSE(UNIX)
  ADD_CUSTOM_COMMAND(
    COMMENT "make package dependencies only implemented in unix"
    TARGET  dep
  )
ENDIF(UNIX)

ADD_CUSTOM_TARGET (check)
IF (UNIX)
  ADD_CUSTOM_COMMAND(
    COMMENT "run tests"
    COMMAND make --directory=../test/ -f Makefile.classic NORMALIZ=${CMAKE_CURRENT_BINARY_DIR}/normaliz
    TARGET  check
    DEPENDS normaliz
  )
  IF(NMZ_COCOA AND COCOA_FOUND)
      ADD_CUSTOM_COMMAND(
      COMMENT "run Int tests"
      COMMAND make --directory=../test/ -f Makefile.classic testInt NORMALIZ=${CMAKE_CURRENT_BINARY_DIR}/normaliz
      TARGET  check
      DEPENDS normaliz
  )
  ENDIF(NMZ_COCOA AND COCOA_FOUND)
ELSE(UNIX)
  ADD_CUSTOM_COMMAND(
    COMMENT "check only implemented in unix"
    TARGET  check
  )
ENDIF(UNIX)
