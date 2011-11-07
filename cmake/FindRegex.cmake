# - Try to find Regex
# Once done this will define
#
#  REGEX_FOUND - system has BZip2
#  REGEX_INCLUDE_DIR - the BZip2 include directory
#  REGEX_LIBRARIES - Link these to use BZip2
#  REGEX_DEFINITIONS - Compiler switches required for using BZip2
#  REGEX_NEED_PREFIX - this is set if the functions are prefixed with BZ2_

# Copyright (c) 2011, Jose Fernandez, jc.fernandez.navarro@gmail.com
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.


IF (REGEX_INCLUDE_DIR AND REGEX_LIBRARIES)
    SET(REGEX_FIND_QUIETLY TRUE)
ENDIF (REGEX_INCLUDE_DIR AND REGEX_LIBRARIES)

FIND_PATH(REGEX_INCLUDE_DIR regex.h )

FIND_LIBRARY(REGEX_LIBRARIES NAMES regex libregex.dll.a libregex.a )

# handle the QUIETLY and REQUIRED arguments and set REGEX_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(REGEX DEFAULT_MSG REGEX_LIBRARIES REGEX_INCLUDE_DIR)


MARK_AS_ADVANCED(REGEX_INCLUDE_DIR REGEX_LIBRARIES)