# - Try to find Glob
# Once done this will define
#
#  GLOB_FOUND - system has Glob
#  GLOB_INCLUDE_DIR - the Glob include directory
#  GLOB_LIBRARIES - Link these to use Glob
#  GLOB_DEFINITIONS - Compiler switches required for using Glob
#  GLOB_NEED_PREFIX - this is set if the functions are prefixed with Glob

# Copyright (c) 2011, Jose Fernandez, jc.fernandez.navarro@gmail.com
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.


IF (GLOB_INCLUDE_DIR AND GLOB_LIBRARIES)
    SET(GLOB_FIND_QUIETLY TRUE)
ENDIF (GLOB_INCLUDE_DIR AND GLOB_LIBRARIES)

FIND_PATH(GLOB_INCLUDE_DIR glob.h )

FIND_LIBRARY(GLOB_LIBRARIES NAMES glob libglob.dll.a libglob.a )

# handle the QUIETLY and REQUIRED arguments and set GLOB_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GLOB DEFAULT_MSG GLOB_LIBRARIES GLOB_INCLUDE_DIR)


MARK_AS_ADVANCED(GLOB_INCLUDE_DIR GLOB_LIBRARIES)