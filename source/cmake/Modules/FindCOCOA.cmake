# - Try to find the COCOA libraries
# This module defines:
#  COCOA_FOUND             - system has COCOA lib
#  COCOA_INCLUDE_DIR       - the COCOA include directory
#  COCOA_LIBRARY_DIR     - directory where the COCOA libraries are located
#  COCOA_LIBRARY         - Link these to use COCOA

include(FindPackageHandleStandardArgs)

if(COCOA_INCLUDE_DIR)
  set(COCOA_in_cache TRUE)
else()
  set(COCOA_in_cache FALSE)
endif()
if(NOT COCOA_LIBRARY)
  set(COCOA_in_cache FALSE)
endif()

# Is it already configured?
if (COCOA_in_cache)

  set(COCOA_FOUND TRUE)

else()

  find_path(COCOA_INCLUDE_DIR
            NAMES CoCoA/library.H
            HINTS ENV COCOA_INC_DIR
                  ENV COCOA_DIR
            PATH_SUFFIXES include
  	        DOC "The directory containing the COCOA header files"
           )

  find_library(COCOA_LIBRARY NAMES lib/libcocoa.a libcocoa
    HINTS ENV COCOA_LIB_DIR
          ENV COCOA_DIR
    PATH_SUFFIXES lib
    DOC "Path to the static COCOA library"
    )

  if ( COCOA_LIBRARY )
    get_filename_component(COCOA_LIBRARY_DIR ${COCOA_LIBRARY} PATH CACHE )
  endif()

  # Attempt to load a user-defined configuration for COCOA if couldn't be found
  if ( NOT COCOA_INCLUDE_DIR OR NOT COCOA_LIBRARY_DIR )
    include( COCOAConfig OPTIONAL )
  endif()

  find_package_handle_standard_args(COCOA "DEFAULT_MSG" COCOA_LIBRARY COCOA_INCLUDE_DIR)

endif()
