# - Try to find Mman
# Once done this will define
#
#  MMAN_FOUND - system has Mman
#  MMAN_INCLUDE_DIR - the Mman include directory
#  MMAN_LIBRARIES - Link these to use Mman
#  MMAN_DEFINITIONS - Compiler switches required for using Mman
#  MMAN_NEED_PREFIX - this is set if the functions are prefixed with Mman

# Copyright (c) 2011, Jose Fernandez, jc.fernandez.navarro@gmail.com
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.


IF (MMAN_INCLUDE_DIR AND MMAN_LIBRARIES)
    SET(MMAN_FIND_QUIETLY TRUE)
ENDIF (MMAN_INCLUDE_DIR AND MMAN_LIBRARIES)

FIND_PATH(MMAN_INCLUDE_DIR mman.h )

FIND_LIBRARY(MMAN_LIBRARIES NAMES mman libmman.dll.a libmman.a )

# handle the QUIETLY and REQUIRED arguments and set MMAN_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MMAN DEFAULT_MSG MMAN_LIBRARIES MMAN_INCLUDE_DIR)


MARK_AS_ADVANCED(MMAN_INCLUDE_DIR MMAN_LIBRARIES)