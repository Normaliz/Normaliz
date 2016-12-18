# - Try to find the GMP libraries for MIC
# This module defines:
#  MIC_GMP_FOUND             - system has GMP lib
#  MIC_GMP_INCLUDE_DIR       - the GMP include directory
#  MIC_GMP_LIBRARIES_DIR     - directory where the GMP libraries are located
#  MIC_GMP_LIBRARIES         - Link these to use GMP
#  MIC_GMP_IN_CGAL_AUXILIARY - TRUE if the GMP found is the one distributed with CGAL in the auxiliary folder

# TODO: support MacOSX

include(FindPackageHandleStandardArgs)
#include(CGAL_GeneratorSpecificSettings)

if(MIC_GMP_INCLUDE_DIR)
  set(MIC_GMP_in_cache TRUE)
else()
  set(MIC_GMP_in_cache FALSE)
endif()
if(NOT MIC_GMP_LIBRARIES)
  set(MIC_GMP_in_cache FALSE)
endif()

# Is it already configured?
if (MIC_GMP_in_cache)

  set(MIC_GMP_FOUND TRUE)

else()

  find_path(MIC_GMP_INCLUDE_DIR
            NAMES gmp.h
            HINTS ENV MIC_GMP_INC_DIR
                  ENV MIC_GMP_DIR
#                  ${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/gmp/include
            PATH_SUFFIXES include
  	        DOC "The directory containing the GMP header files"
           )

  if ( MIC_GMP_INCLUDE_DIR STREQUAL "${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/gmp/include" )
    cache_set( MIC_GMP_IN_CGAL_AUXILIARY TRUE )
  endif()

  find_library(MIC_GMP_LIBRARIES NAMES gmp libgmp-10
    HINTS ENV MIC_GMP_LIB_DIR
          ENV MIC_GMP_DIR
#          ${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/gmp/lib
    PATH_SUFFIXES lib
    DOC "Path to the GMP library"
    )
  find_library(MIC_GMP_STATIC_LIBRARIES NAMES libgmp.a
    HINTS ENV MIC_GMP_LIB_DIR
          ENV MIC_GMP_DIR
    PATH_SUFFIXES lib
    DOC "Path to the static GMP library"
    )

  if ( MIC_GMP_LIBRARIES )
    get_filename_component(MIC_GMP_LIBRARIES_DIR ${MIC_GMP_LIBRARIES} PATH CACHE )
  endif()

  # Attempt to load a user-defined configuration for GMP if couldn't be found
  if ( NOT MIC_GMP_INCLUDE_DIR OR NOT MIC_GMP_LIBRARIES_DIR )
    include( MIC_GMPConfig OPTIONAL )
  endif()

  find_package_handle_standard_args(MIC_GMP "DEFAULT_MSG" MIC_GMP_LIBRARIES MIC_GMP_INCLUDE_DIR)

endif()
