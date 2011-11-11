# - Try to find PSapi
# Once done this will define
#
#  PSAPI_FOUND - system has PSapi
#  PSAPI_INCLUDE_DIR - the Mman include directory
#  PSAPI_LIBRARIES - Link these to use PSapi
#  PSAPI_DEFINITIONS - Compiler switches required for using PSapi
#  PSAPI_NEED_PREFIX - this is set if the functions are prefixed with PSapi

# Copyright (c) 2011, Jose Fernandez, jc.fernandez.navarro@gmail.com
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.


IF (PSAPI_INCLUDE_DIR AND PSAPI_LIBRARIES)
    SET(PSAPI_FIND_QUIETLY TRUE)
ENDIF (PSAPI_INCLUDE_DIR AND PSAPI_LIBRARIES)

FIND_PATH(PSAPI_INCLUDE_DIR psapi.h )

FIND_LIBRARY(PSAPI_LIBRARIES NAMES psapi libpsapi.dll.a libpsapi.a )

# handle the QUIETLY and REQUIRED arguments and set PSAPI_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PSAPI DEFAULT_MSG PSAPI_LIBRARIES PSAPI_INCLUDE_DIR)


MARK_AS_ADVANCED(PSAPI_INCLUDE_DIR PSAPI_LIBRARIES)