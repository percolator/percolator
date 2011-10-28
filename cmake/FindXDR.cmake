# Try to find the XDR library and set some required variables
#
# Once run this will define:
#
# XDR_FOUND                                     = system has XDR lib
#
# XDR_LIBRARIES                 = full path to the libraries, if required
#

INCLUDE(CheckIncludeFileCXX)
INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckLibraryExists)

## First try to find the required header files (rpc/types.h, rpc/xdr.h)

#check for the XDR functions: their interface and the libraries they're hidden in.
CHECK_INCLUDE_FILE_CXX(rpc/types.h XDR_HAS_RPC_TYPES_H)
IF(NOT XDR_HAS_RPC_TYPES_H)
    MESSAGE(SEND_ERROR "Cannot find RPC headers (rpc/types.h)")
ELSE()
    CHECK_CXX_SOURCE_COMPILES("#include <rpc/types.h>
			       #include <rpc/xdr.h>

			      int main(int /*argc*/, char** /*argv*/)
			      {
				return 0;
			      }" XDR_HAS_RPC_XDR_H)
    IF(NOT XDR_HAS_RPC_XDR_H)
      MESSAGE(SEND_ERROR "Cannot find RPC headers (rpc/xdr.h)")
    ENDIF()
ENDIF()
	
IF (XDR_HAS_RPC_TYPES_H AND XDR_HAS_RPC_XDR_H)
     ## Now let's see if we need an extra lib to compile it

  find_library(XDR_LIBRARIES NAMES libportablexdr.a libportablexdr.dll.a libportablexdr.la HINTS
     /usr/i686-pc-mingw32/sys-root/mingw/lib/
     /usr/i586-mingw32msvc/sys-root/mingw/lib/
     /usr/i586-mingw32msvc/lib/
     /usr/i686-pc-mingw32/lib/
     C:/MinGW/lib/
     /mingw/lib   
  )

  message(STATUS "XDR_LIBRARIES=${XDR_LIBRARIES}")

  SET(XDR_INT_FOUND FALSE)
  
  #check functions exist
  CHECK_LIBRARY_EXISTS(${XDR_LIBRARIES} xdr_u_int "" XDR_INT_FOUND)
  CHECK_LIBRARY_EXISTS(${XDR_LIBRARIES} xdr_bool "" XDR_INT_FOUND)
  CHECK_LIBRARY_EXISTS(${XDR_LIBRARIES} xdr_double "" XDR_INT_FOUND)
  CHECK_LIBRARY_EXISTS(${XDR_LIBRARIES} xdr_uint32_t "" XDR_INT_FOUND)
  CHECK_LIBRARY_EXISTS(${XDR_LIBRARIES} xdr_int32_t "" XDR_INT_FOUND)
  CHECK_LIBRARY_EXISTS(${XDR_LIBRARIES} xdr_opaque "" XDR_INT_FOUND)
  
  IF(NOT XDR_INT_FOUND OR NOT XDR_LIBRARIES)
    MESSAGE(SEND_ERROR "Could not locate xdr symbols")
  ELSE()
    SET(XDR_FOUND TRUE)
  ENDIF()

ENDIF()
