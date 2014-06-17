#                                               -*- cmake -*-
#
#  UseLibnormaliz.cmake
#
#  Copyright (C) 2014 Christof Soeger <csoeger@uos.de>
#

#FIND_PACKAGE(LibFTDI1 NO_MODULE REQUIRED)
#INCLUDE(${LIBFTDI_USE_FILE})

add_definitions     ( ${LIBNORAMLIZ_DEFINITIONS} )
include_directories ( ${LIBNORAMLIZ_INCLUDE_DIRS} )
link_directories    ( ${LIBNORAMLIZ_LIBRARY_DIRS} )
