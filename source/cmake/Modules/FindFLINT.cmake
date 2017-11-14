# - Try to find the FLINT libraries
# This module defines:
#  FLINT_FOUND             - system has FLINT lib
#  FLINT_INCLUDE_DIR       - the FLINT include directory
#  FLINT_LIBRARY_DIR     - directory where the FLINT libraries are located
#  FLINT_LIBRARY         - Link these to use FLINT

include(FindPackageHandleStandardArgs)

if(FLINT_INCLUDE_DIR)
  set(FLINT_in_cache TRUE)
else()
  set(FLINT_in_cache FALSE)
endif()
if(NOT FLINT_LIBRARY)
  set(FLINT_in_cache FALSE)
endif()

# Is it already configured?
if (FLINT_in_cache)

  set(FLINT_FOUND TRUE)

else()

  find_path(FLINT_INCLUDE_DIR
            NAMES  flint/flint.h
            HINTS ENV FLINT_INC_DIR
                  ENV FLINT_DIR
            PATH_SUFFIXES include
  	        DOC "The directory containing the FLINT header files"
           )

  find_library(FLINT_LIBRARY NAMES lib/libflint.so libflint.so
    HINTS ENV FLINT_LIB_DIR
          ENV FLINT_DIR
    PATH_SUFFIXES lib
    DOC "Path to the static FLINT library"
    )

  if ( FLINT_LIBRARY )
    get_filename_component(FLINT_LIBRARY_DIR ${FLINT_LIBRARY} PATH CACHE )
  endif()

  # Attempt to load a user-defined configuration for FLINT if couldn't be found
  if ( NOT FLINT_INCLUDE_DIR OR NOT FLINT_LIBRARY_DIR )
    include( FLINTConfig OPTIONAL )
  endif()

  find_package_handle_standard_args(FLINT "DEFAULT_MSG" FLINT_LIBRARY FLINT_INCLUDE_DIR)

endif()
