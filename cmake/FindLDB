# - Try to find LevelDB
# Once done this will define
#
# LDB_FOUND - system has LevelDB
# LDB_INCLUDE_DIR - the LevelDB include directory
# LDB_LIBRARIES - Link these to use LevelDB
# LDB_DEFINITIONS - Compiler switches required for using LevelDB
# LDB_NEED_PREFIX - this is set if the functions are prefixed with LevelDB

# Copyright (c) 2011, Jose Fernandez, jc.fernandez.navarro@gmail.com
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.


IF (LDB_INCLUDE_DIR AND LDB_LIBRARIES)
    SET(LDB_FIND_QUIETLY TRUE)
ENDIF (LDB_INCLUDE_DIR AND LDB_LIBRARIES)

FIND_PATH(LDB_INCLUDE_DIR NAMES db.h PATH_SUFFIXES leveldb HINTS
   /usr/include
   /usr/local/include
   $ENV{LEVELDB}
   $ENV{LEVELDB}/include
   )

FIND_LIBRARY(LDB_LIBRARIES NAMES leveldb leveldb.dll.a leveldb.a HINTS
   /usr/lib
   /usr/local/lib
   $ENV{LEVELDB}
   $ENV{LEVELDB}/lib )

# handle the QUIETLY and REQUIRED arguments and set BZip2_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LDB DEFAULT_MSG LDB_LIBRARIES LDB_INCLUDE_DIR)

MARK_AS_ADVANCED(LDB_INCLUDE_DIR LDB_LIBRARIES)
