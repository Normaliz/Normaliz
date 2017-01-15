# libQnormaliz source files
SET(sources Qcone_property.cpp libQnormaliz-templated.cpp Qcone_helper.cpp)
# libnormaliz header files
SET(headers Qcone.h Qcone_property.h Qconvert.h Qgeneral.h Qinteger.h libQnormaliz.h Qmap_operations.h Qmatrix.h Qmy_omp.h Qnormaliz_exception.h Qsublattice_representation.h Qvector_operations.h Qversion.h)

# add the necessary include directories
INCLUDE_DIRECTORIES(..)
FIND_PACKAGE(Boost REQUIRED)
#IF (Boost_FOUND)
#  MESSAGE(STATUS "Boost found with include dir: ${Boost_INCLUDE_DIR}")
#ENDIF()
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

SET (QNMZ_BUILD_BOTH "NO" CACHE BOOL "build shared AND static libQnormaliz")

IF (QNMZ_BUILD_BOTH OR NOT QNMZ_BUILD_STATIC)
  # create the shared library
  ADD_LIBRARY(Qnormaliz-shared SHARED ${sources})
  SET_TARGET_PROPERTIES(Qnormaliz-shared PROPERTIES OUTPUT_NAME Qnormaliz VERSION ${VERSION} SOVERSION 0)
  # set the dependency on gmp
  TARGET_LINK_LIBRARIES (Qnormaliz-shared ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
  INSTALL(TARGETS Qnormaliz-shared
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib${LIB_SUFFIX}
    RUNTIME DESTINATION bin
    COMPONENT sharedlibs
  )
ENDIF ()

IF (QNMZ_BUILD_BOTH OR QNMZ_BUILD_STATIC)
  #create the static library
  ADD_LIBRARY(Qnormaliz-static STATIC ${sources})
  SET_TARGET_PROPERTIES(Qnormaliz-static PROPERTIES OUTPUT_NAME Qnormaliz VERSION ${VERSION} SOVERSION 0)
  INSTALL(TARGETS Qnormaliz-static ARCHIVE DESTINATION lib${LIB_SUFFIX})
ENDIF ()

INSTALL(FILES ${headers} DESTINATION include/${PROJECT_NAME} COMPONENT headers)

SET(LIBQNORMALIZ_INCLUDE_DIR ${includedir})
SET(LIBQNORMALIZ_INCLUDE_DIRS ${LIBQNORMALIZ_INCLUDE_DIR})
SET(LIBQNORMALIZ_ROOT_DIR ${prefix})
SET(LIBQNORMALIZ_VERSION_STRING ${VERSION_STRING})
SET(LIBQNORMALIZ_VERSION_MAJOR ${VERSION_MAJOR})
SET(LIBQNORMALIZ_VERSION_MINOR ${VERSION_MINOR})
SET(LIBQNORMALIZ_VERSION_PATCH ${VERSION_PATCH})
CONFIGURE_FILE(Qversion.h.in ${CMAKE_CURRENT_SOURCE_DIR}/Qversion.h @ONLY)
SET(LIBQNORMALIZ_USE_FILE ${CMAKE_INSTALL_PREFIX}/${LIBQNORMALIZ_CMAKE_CONFIG_DIR}/UseLibQnormaliz.cmake)
CONFIGURE_FILE(../LibQnormalizConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/LibQnormalizConfig.cmake @ONLY)
#INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/LibQnormalizConfig.cmake ../UseLibnormaliz.cmake DESTINATION ${LIBQNORMALIZ_CMAKE_CONFIG_DIR})

#ADD_SUBDIRECTORY(examples)
