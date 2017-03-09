# - Try to find the NAUTY libraries
# This module defines:
#  NAUTY_FOUND             - system has NAUTY lib
#  NAUTY_INCLUDE_DIR       - the NAUTY include directory
#  NAUTY_LIBRARY_DIR     - directory where the NAUTY libraries are located
#  NAUTY_LIBRARY         - Link these to use NAUTY

include(FindPackageHandleStandardArgs)

if(NAUTY_INCLUDE_DIR)
  set(NAUTY_in_cache TRUE)
else()
  set(NAUTY_in_cache FALSE)
endif()
if(NOT NAUTY_LIBRARY)
  set(NAUTY_in_cache FALSE)
endif()

# Is it already configured?
if (NAUTY_in_cache)

  set(NAUTY_FOUND TRUE)

else()

  find_path(NAUTY_INCLUDE_DIR
            NAMES nauty/nauty.h
            HINTS ENV NAUTY_INC_DIR
                  ENV NAUTY_DIR
            PATH_SUFFIXES include
  	        DOC "The directory containing the NAUTY header files"
           )

  find_library(NAUTY_LIBRARY NAMES nauty
    HINTS ENV NAUTY_LIB_DIR
          ENV NAUTY_DIR
    PATH_SUFFIXES lib
    DOC "Path to the NAUTY library"
    )
  find_library(NAUTY_STATIC_LIBRARY NAMES libnauty.a
    HINTS ENV NAUTY_LIB_DIR
          ENV NAUTY_DIR
    PATH_SUFFIXES lib
    DOC "Path to the static NAUTY library"
    )

  if ( NAUTY_LIBRARY )
    get_filename_component(NAUTY_LIBRARY_DIR ${NAUTY_LIBRARY} PATH CACHE )
  endif()

  # Attempt to load a user-defined configuration for NAUTY if couldn't be found
  if ( NOT NAUTY_INCLUDE_DIR OR NOT NAUTY_LIBRARY_DIR )
    include( NAUTYConfig OPTIONAL )
  endif()

  find_package_handle_standard_args(NAUTY "DEFAULT_MSG" NAUTY_LIBRARY NAUTY_INCLUDE_DIR)

endif()
