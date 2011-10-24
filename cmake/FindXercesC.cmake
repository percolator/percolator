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
   /usr/include
   /usr/local/include
   $ENV{XERCESCROOT}/src
   $ENV{XERCESCROOT}/include
   $ENV{CODESYNTH}/include
   $ENV{XSDDIR}/include
   ${XERCES_DIR}/include
   $ENV{XERCES_DIR}/include
   ${XERCESCROOT}/include
   ${CODESYNTH}/include
   ${XSDDIR}/include
   PATH_SUFFIXES xerces-c xercesc include
   )

FIND_LIBRARY(XERCESC_LIBRARIES NAMES xerces-c xerces-c_3 xerces-c_2 xerces-c_static xerces-c_static_3 xerces-c_static_2 libxerces-c
   HINTS
   "[HKEY_CURRENT_USER\\software\\xerces-c\\lib]"
   "[HKEY_CURRENT_USER\\xerces-c\\lib]"
   /usr/lib
   /usr/local/lib
   $ENV{XERCESCROOT}/lib
   $ENV{CODESYNTH}/lib
   $ENV{XSDDIR}/lib
   ${XERCES_DIR}/lib
   $ENV{XERCES_DIR}/lib
   ${XERCESCROOT}/lib
   ${CODESYNTH}/lib
   ${XSDDIR}/lib
   PATH_SUFFIXES lib64 lib32 lib
   )

if( MSVC )
# make sure we use the right compiler version

	if( MSVC80 )
		set( XERCESC_COMPILER_PREFIX vc-8.0 )
	elseif( MSVC90 )
		set( XERCESC_COMPILER_PREFIX vc-9.0 )
	elseif( MSVC10 )
		set( XERCESC_COMPILER_PREFIX vc-10.0 )
	else( MSVC10 )
		message( FATAL_ERROR "Unknown MSVC version!" )
	endif( MSVC80 ) 
		
endif( MSVC )

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