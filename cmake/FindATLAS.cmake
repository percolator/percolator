# - Try to find atlas
# Once done this will define
#
# This module defines
# XSD_INCLUDE_DIR, where to find elements.hxx, etc.
# XSD_EXECUTABLE, the exe file
# XSD_FOUND, If false, don't try to use it.

FIND_PATH(ATLAS_INCLUDE_DIR atlas/clapack.h
PATHS
  $ENV{ATLASDIR}
  /usr/local
  /usr
PATH_SUFFIXES
  include
  libatlas
  atlas
  NO_CMAKE_FIND_ROOT_PATH
)


# if the include and the program are found then we have it
IF(ATLAS_INCLUDE_DIR)
  SET( ATLAS_FOUND "YES" )
ENDIF(ATLAS_INCLUDE_DIR)


MARK_AS_ADVANCED(
  ATLAS_INCLUDE_DIR
) 
