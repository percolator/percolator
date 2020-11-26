# - Try to find XercesC
# Once done this will define
#
# This module defines
# XERCESC_INCLUDE_DIR, where to find xercesc.h, etc.
# XERCESC_LIBRARIES, the libraries to link against to use xercesc.
# XERCESC_FOUND, If false, don't try to use xercesc.
# XERCESC_VERSION, version of the library found

FIND_PATH(XERCESC_INCLUDE_DIR xercesc/dom/DOM.hpp xercesc/parsers/SAXParser.hpp
   HINTS
   "[HKEY_CURRENT_USER\\software\\xerces-c\\src]"
   "[HKEY_CURRENT_USER\\xerces-c\\src]"
   $ENV{XERCESCROOT}/src
   $ENV{XERCESCROOT}
   /usr/include
   /usr/local/include
   "[HKEY_CURRENT_USER\\software\\xsd\\include]"
   "[HKEY_CURRENT_USER]\\xsd\\include]"
   PATH_SUFFIXES xerces-c xercesc include
   )

FIND_LIBRARY(XERCESC_LIBRARIES NAMES xerces-c xerces-c_3 xerces-c_static_3 libxerces-c
   HINTS
   "[HKEY_CURRENT_USER\\software\\xerces-c\\lib]"
   "[HKEY_CURRENT_USER\\xerces-c\\lib]"
   $ENV{XERCESCROOT}/
   $ENV{XERCESCROOT}/lib
   $ENV{XERCESCROOT}/src
   /usr/
   /usr/lib
   /usr/local/lib
   PATH_SUFFIXES src lib64 lib32 lib
   )

if( MSVC )
# make sure we use the right compiler version
	if( MSVC80 )
		set( XERCESC_COMPILER_PREFIX vc-8.0 )
	elseif( MSVC90 )
		set( XERCESC_COMPILER_PREFIX vc-9.0 )
	elseif( MSVC10 )
		set( XERCESC_COMPILER_PREFIX vc-10.0 )  
	elseif( MSVC11 )
		set( XERCESC_COMPILER_PREFIX vc-10.0 )
		message( WARNING "No Xerces available for VS2012, using the version for VS2010 instead")
	elseif( MSVC12 )
		set( XERCESC_COMPILER_PREFIX vc-10.0 )
		message( WARNING "No Xerces available for VS2013, using the version for VS2010 instead")
  elseif( MSVC14 )
		set( XERCESC_COMPILER_PREFIX vc-10.0 )
		message( WARNING "No Xerces available for VS2013, using the version for VS2010 instead")
	else( MSVC10 )
		message( FATAL_ERROR "Unknown MSVC version!" )
	endif( MSVC80 ) 
		
endif( MSVC )

# if the include a the library are found then we have it
## CHECK  FOR CORRECT VERSION OF XERCESC AND THE FUNCTIONS INCLUDED
SET( XERCESC_FOUND 0 )
IF(XERCESC_INCLUDE_DIR)
  IF(XERCESC_LIBRARIES)
    SET( XERCESC_FOUND 1 )
    FIND_PATH(XERCESC_XVERHPPPATH NAMES XercesVersion.hpp PATHS ${XERCESC_INCLUDE_DIR} PATH_SUFFIXES xercesc/util)
    IF ( ${XERCESC_XVERHPPPATH} STREQUAL XERCESC_XVERHPPPATH-NOTFOUND )
      SET(XERCES_VERSION "0")
    ELSE( ${XERCESC_XVERHPPPATH} STREQUAL XERCESC_XVERHPPPATH-NOTFOUND )
      FILE(READ ${XERCESC_XVERHPPPATH}/XercesVersion.hpp XVERHPP)
      STRING(REGEX MATCHALL "\n *#define XERCES_VERSION_MAJOR +[0-9]+" XVERMAJ ${XVERHPP})
      STRING(REGEX MATCH "\n *#define XERCES_VERSION_MINOR +[0-9]+" XVERMIN ${XVERHPP})
      STRING(REGEX MATCH "\n *#define XERCES_VERSION_REVISION +[0-9]+" XVERREV ${XVERHPP})
      STRING(REGEX REPLACE "\n *#define XERCES_VERSION_MAJOR +" "" XVERMAJ ${XVERMAJ})
      STRING(REGEX REPLACE "\n *#define XERCES_VERSION_MINOR +" "" XVERMIN ${XVERMIN})
      STRING(REGEX REPLACE "\n *#define XERCES_VERSION_REVISION +" "" XVERREV ${XVERREV})
      SET(XERCESC_VERSION ${XVERMAJ}.${XVERMIN}.${XVERREV})
      if(${XVERMAJ} LESS 3)
	message(FATAL_ERROR "The version of Xerces-c found : " ${XERCESC_VERSION} " is too old ")
      endif()
    ENDIF ( ${XERCESC_XVERHPPPATH} STREQUAL XERCESC_XVERHPPPATH-NOTFOUND )   
  ENDIF(XERCESC_LIBRARIES)
ENDIF(XERCESC_INCLUDE_DIR)



# include(CheckCXXSourceCompiles)
# set(CMAKE_REQUIRED_INCLUDES ${XERCESC_INCLUDE_DIR})
# set(CMAKE_REQUIRED_LIBRARIES ${XERCESC_LIBRARIES})
# CHECK_CXX_SOURCE_COMPILES("    #include <xercesc/dom/DOM.hpp>
# 				 #include <xercesc/util/XMLString.hpp>
# 				 #include <xercesc/parsers/XercesDOMParser.hpp>
# 				 #include <xercesc/sax/HandlerBase.hpp>
# 				 #include <xercesc/util/PlatformUtils.hpp>
#                                using namespace xercesc;
# 				 int main(int /*argc*/, char** /*argv*/)
# 				 {
# 			           xercesc::XMLPlatformUtils::Initialize();
# 				   return 0;
# 			         }" 
# 			        HAS_XERCESC_XDR_C)
#     IF(NOT HAS_XERCESC_XDR_C)
#       MESSAGE(FATAL_ERROR "Cannot compile Xerces-c code")
#     ENDIF()
# ENDIF()

MARK_AS_ADVANCED(
  XERCESC_INCLUDE_DIR
  XERCESC_LIBRARIES
  XERCESC_VERSION
)
