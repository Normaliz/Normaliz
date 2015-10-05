# - Try to find the scipoptsuite library
# This module defines:
#  SCIP_FOUND          - system has SCIP lib
#  SCIP_INCLUDE_DIR    - the SCIP include directory
#  SCIP_LIBRARY_DIR    - directory where the SCIP libraries are located
#  SCIP_LIBRARIES      - Link these to use SCIP, it contains
#  SCIP_LIBRARY        - libscipopt
#  Readline_LIBRARY    - libreadline
#  Z_LIBRARY           - libz

#TODO version is fixed!

INCLUDE(FindPackageHandleStandardArgs)

IF(SCIP_INCLUDE_DIR AND SCIP_LIBRARY)
  SET(SCIP_in_cache TRUE)
ELSE()
  SET(SCIP_in_cache FALSE)
ENDIF()
IF(NOT SCIP_LIBRARY)
  SET(SCIP_in_cache FALSE)
ENDIF()

# Is it already configured?
IF (SCIP_in_cache)

  SET(SCIP_FOUND TRUE)

ELSE()

  FIND_PATH(SCIP_INCLUDE_DIR
            NAMES scip/scip.h
            HINTS ENV SCIP_INC_DIR
                  ENV SCIP_DIR
            PATH_SUFFIXES include src scip-3.1.1/src scip-3.2.0/src
  	        DOC "The directory containing the SCIP header files"
  )

  FIND_LIBRARY(SCIP_LIBRARY 
    NAMES libscipopt-3.2.0.linux.x86_64.gnu.opt.a libscipopt-3.1.1.linux.x86_64.gnu.opt.a libscipopt-3.2.0.darwin.x86_64.gnu.opt.a  libscipopt-3.1.1.darwin.x86_64.gnu.opt.a
    HINTS ENV SCIP_LIB_DIR
          ENV SCIP_DIR
    PATH_SUFFIXES lib
    DOC "Path to the SCIP library"
  )

  IF ( SCIP_LIBRARY )
    GET_FILENAME_COMPONENT( SCIP_LIBRARY_DIR ${SCIP_LIBRARY} PATH CACHE )
  ENDIF()
  #SET(SCIP_LIBRARIES "${SCIP_LIBRARY}")
  SET(SCIP_LIBRARIES "${SCIP_LIBRARY} ${Readline_LIBRARY} ${Z_LIBRARY}")

  FIND_PACKAGE_HANDLE_STANDARD_ARGS(SCIP DEFAULT_MSG SCIP_INCLUDE_DIR SCIP_LIBRARY)
ENDIF()
