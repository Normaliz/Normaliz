# - Try to find the MPFR libraries
# This module defines:
#  MPFR_FOUND             - system has MPFR lib
#  MPFR_INCLUDE_DIR       - the MPFR include directory
#  MPFR_LIBRARY_DIR     - directory where the MPFR libraries are located
#  MPFR_LIBRARY         - Link these to use MPFR

include(FindPackageHandleStandardArgs)

if(MPFR_INCLUDE_DIR)
  set(MPFR_in_cache TRUE)
else()
  set(MPFR_in_cache FALSE)
endif()
if(NOT MPFR_LIBRARY)
  set(MPFR_in_cache FALSE)
endif()

# Is it already configured?
if (MPFR_in_cache)

  set(MPFR_FOUND TRUE)

else()

  find_path(MPFR_INCLUDE_DIR
            NAMES mpfr.h
            HINTS ENV MPFR_INC_DIR
                  ENV MPFR_DIR
            PATH_SUFFIXES include
  	        DOC "The directory containing the MPFR header files"
           )

  find_library(MPFR_LIBRARY NAMES lib/libmpfr.so libmpfr.so
    HINTS ENV MPFR_LIB_DIR
          ENV MPFR_DIR
    PATH_SUFFIXES lib
    DOC "Path to the shared MPFR library"
    )
    
    find_library(MPFR_STATIC_LIBRARY NAMES b/libmpfr.a libmpfr.a
    HINTS ENV FLINT_LIB_DIR
          ENV FLINT_LIB_DIR
    PATH_SUFFIXES lib
    DOC "Path to the static MPFR library"
    )

  if ( MPFR_LIBRARY )
    get_filename_component(MPFR_LIBRARY_DIR ${MPFR_LIBRARY} PATH CACHE )
  endif()

  # Attempt to load a user-defined configuration for MPFR if couldn't be found
  if ( NOT MPFR_INCLUDE_DIR OR NOT MPFR_LIBRARY_DIR )
    include( MPFRConfig OPTIONAL )
  endif()

  find_package_handle_standard_args(MPFR "DEFAULT_MSG" MPFR_LIBRARY MPFR_INCLUDE_DIR)

endif()
