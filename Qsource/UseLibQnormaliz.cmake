#                                               -*- cmake -*-
#
#  UseQLibnormaliz.cmake
#
#  Copyright (C) 2014 Christof Soeger <csoeger@uos.de>
#

#FIND_PACKAGE(LibFTDI1 NO_MODULE REQUIRED)
#INCLUDE(${LIBFTDI_USE_FILE})

add_definitions     ( ${LIBQNORAMLIZ_DEFINITIONS} )
include_directories ( ${LIBQNORAMLIZ_INCLUDE_DIRS} )
link_directories    ( ${LIBQNORAMLIZ_LIBRARY_DIRS} )
