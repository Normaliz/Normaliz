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

# The project name and the type of project
PROJECT(libnormaliz)
SET(MAJOR_VERSION 2)
SET(MINOR_VERSION 11)
SET(PATCH_VERSION 1)
SET(VERSION_STRING ${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION})
SET(VERSION ${VERSION_STRING})

# automatically set lib suffix
IF(UNIX AND NOT APPLE AND NOT CMAKE_CROSSCOMPILING AND NOT EXISTS "/etc/debian_version")
  IF(CMAKE_SIZEOF_VOID_P EQUAL 8 AND NOT LIB_SUFFIX)
    SET(LIB_SUFFIX 64)
  ENDIF()
ENDIF()

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

# compiler flags
SET (CXX_WARNING_FLAGS "-Wall -pedantic -Wno-long-long")
SET (CXX_OPTIMIZATION_FLAGS "-O3 -funroll-loops")
SET (NMZ_OPENMP "YES" CACHE BOOL "use openmp")
IF (NMZ_OPENMP)
  FIND_PACKAGE(OpenMP)
  IF (OPENMP_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  ELSE (OPENMP_FOUND)
    MESSAGE(STATUS "Not using OpenMP, compiler does not support it")
    SET (NMZ_OPENMP "NO")
  ENDIF(OPENMP_FOUND)
ENDIF (NMZ_OPENMP)
IF (NOT NMZ_OPENMP)
  SET (CXX_WARNING_FLAGS "${CXX_WARNING_FLAGS} -Wno-unknown-pragmas")
ENDIF (NOT NMZ_OPENMP)

# add compiler flags
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_WARNING_FLAGS} ${CXX_OPTIMIZATION_FLAGS}")

# handle static compilation
SET (NMZ_BUILD_STATIC "NO" CACHE BOOL "build normaliz statically")
IF (NMZ_BUILD_STATIC)
  IF (${APPLE})
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
  ELSEIF (${UNIX})
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  ELSE ()
    MESSAGE (FATAL_ERROR "Static linking not (yet) supported on this system")
  ENDIF ()
ENDIF (NMZ_BUILD_STATIC)

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

ADD_SUBDIRECTORY(libnormaliz)

# add the normaliz executable
ADD_EXECUTABLE (normaliz Normaliz.cpp)

IF (NMZ_BUILD_STATIC)
  TARGET_LINK_LIBRARIES (normaliz normaliz-static)
  TARGET_LINK_LIBRARIES (normaliz ${GMPXX_STATIC_LIBRARIES} ${GMP_STATIC_LIBRARIES})
ELSE (NMZ_BUILD_STATIC)
  TARGET_LINK_LIBRARIES (normaliz normaliz-shared)
  TARGET_LINK_LIBRARIES (normaliz ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
ENDIF (NMZ_BUILD_STATIC)



INSTALL(TARGETS normaliz DESTINATION bin)



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
