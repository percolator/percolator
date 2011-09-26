# - Try to find XercesC
# Once done this will define
#
#  XERCESC_FOUND - System has XercesC
#  XERCESC_INCLUDE_DIR - The XercesC include directory
#  XERCESC_LIBRARY_DIR - The XercesC library dir
#  XERCESC_LIBRARIES - The libraries needed to use XercesC
#  XERCESC_DEFINITIONS - Compiler switches required for using XercesC

# Copyright (c) 2009, Helio Chissini de Castro, <helio@kde.org>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.


IF (XERCESC_INCLUDE_DIR AND XERCESC_LIBRARIES)
   # in cache already
   SET(XercesC_FIND_QUIETLY TRUE)
ENDIF (XERCESC_INCLUDE_DIR AND XERCESC_LIBRARIES)

IF (NOT WIN32)
   # use pkg-config to get the directories and then use these values
   # in the FIND_PATH() and FIND_LIBRARY() calls
   FIND_PACKAGE(PkgConfig)
   PKG_CHECK_MODULES(PC_XERCESC xerces-c)
   SET(XERCESC_DEFINITIONS ${PC_XERCESC_CFLAGS_OTHER})
   SET(XERCESC_LIBRARY_DIR ${PC_XERCESC_LIBRARY_DIRS})
ENDIF (NOT WIN32)

FIND_PATH(XERCESC_INCLUDE_DIR xercesc/dom/DOM.hpp
   HINTS
   ${PC_XERCESC_INCLUDEDIR}
   ${PC_XERCESC_INCLUDE_DIRS}
   /usr/include
   /usr/local/include
   PATH_SUFFIXES xerces-c
   )

FIND_LIBRARY(XERCESC_LIBRARIES NAMES xerces-c xerces-c_3 xerces-c_2 xerces-c_static xerces-c_static_3 xerces-c_static_2 libxerces-c
   HINTS
   ${PC_XERCESC_LIBDIR}
   ${PC_XERCESC_LIBRARY_DIRS}
   /usr/lib
   /usr/local/lib
   )

INCLUDE(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set XERCESC_FOUND to TRUE if 
# all listed variables are TRUE
FIND_PACKAGE_HANDLE_STANDARD_ARGS(XercesC DEFAULT_MSG XERCESC_LIBRARIES XERCESC_INCLUDE_DIR)

MARK_AS_ADVANCED(XERCESC_INCLUDE_DIR XERCESC_LIBRARIES XERCESC_LIBRARY_DIR)
