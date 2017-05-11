# - Try to find the GMPXX libraries for MIC
# This module defines:
#   MIC_GMPXX_FOUND        - system has GMPXX lib
#   MIC_GMPXX_INCLUDE_DIR  - the GMPXX include directory
#   MIC_GMPXX_LIBRARIES    - Libraries needed to use GMPXX

# TODO: support Windows and MacOSX

# GMPXX needs GMP

find_package( MIC_GMP QUIET )

if(MIC_GMP_FOUND)

  if (MIC_GMPXX_INCLUDE_DIR AND MIC_GMPXX_LIBRARIES)
    # Already in cache, be silent
    set(MIC_GMPXX_FIND_QUIETLY TRUE)
  endif()

  find_path(MIC_GMPXX_INCLUDE_DIR NAMES gmpxx.h
            HINTS ENV MIC_GMPXX_INC_DIR
                  ENV MIC_GMPXX_DIR
                  ENV MIC_GMP_INC_DIR
                  ENV MIC_GMP_DIR
                  ${MIC_GMP_INCLUDE_DIR_SEARCH}
            PATH_SUFFIXES include
            DOC "The directory containing the GMPXX include files"
           )

  find_library(MIC_GMPXX_LIBRARIES NAMES gmpxx
               HINTS ENV MIC_GMPXX_LIB_DIR
                     ENV MIC_GMPXX_DIR
                     ENV MIC_GMP_LIB_DIR
                     ENV MIC_GMP_DIR
                     ${MIC_GMP_LIBRARIES_DIR_SEARCH}
               PATH_SUFFIXES lib
               DOC "Path to the GMPXX library"
               )

  find_library(MIC_GMPXX_STATIC_LIBRARIES NAMES libgmpxx.a
               HINTS ENV MIC_GMPXX_LIB_DIR
                     ENV MIC_GMPXX_DIR
                     ENV MIC_GMP_LIB_DIR
                     ENV MIC_GMP_DIR
                     ${MIC_GMP_LIBRARIES_DIR_SEARCH}
               PATH_SUFFIXES lib
               DOC "Path to the static GMPXX library"
               )

  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args(MIC_GMPXX "DEFAULT_MSG" MIC_GMPXX_LIBRARIES MIC_GMPXX_INCLUDE_DIR )

else()

  message( FATAL_ERROR "MIC_GMPXX needs MIC_GMP")

endif()
