# - Try to find XercesC
# Once done this will define
#
# This module defines
# XERCESC_INCLUDE_DIR, where to find ptlib.h, etc.
# XERCESC_LIBRARIES, the libraries to link against to use pwlib.
# XERCESC_FOUND, If false, don't try to use pwlib.

FIND_PATH(XERCESC_INCLUDE_DIR xercesc/dom/DOM.hpp xercesc/parsers/SAXParser.hpp
   HINTS
   "[HKEY_CURRENT_USER\\software\\xerces-c\\src]"
   "[HKEY_CURRENT_USER\\xerces-c\\src]"
   $ENV{XERCESCROOT}/src/
   $ENV{XERCESCROOT}/include/
   /usr/include
   /usr/local/include
   PATH_SUFFIXES xerces-c include
   )

FIND_LIBRARY(XERCESC_LIBRARIES NAMES xerces-c xerces-c_3 xerces-c_2 xerces-c_static xerces-c_static_3 xerces-c_static_2 libxerces-c
   HINTS
   "[HKEY_CURRENT_USER\\software\\xerces-c\\lib]"
   "[HKEY_CURRENT_USER\\xerces-c\\lib]"
   $ENV{XERCESCROOT}/lib
   /usr/lib
   /usr/local/lib
   PATH_SUFFIXES lib64 lib32 lib
   )

# if the include a the library are found then we have it
SET( XERCESC_FOUND 0 )
IF(XERCESC_INCLUDE_DIR)
  IF(XERCESC_LIBRARIES)
    SET( XERCESC_FOUND 1 )
  ENDIF(XERCESC_LIBRARIES)
ENDIF(XERCESC_INCLUDE_DIR)

MARK_AS_ADVANCED(
  XERCESC_INCLUDE_DIR
  XERCESC_LIBRARIES
)