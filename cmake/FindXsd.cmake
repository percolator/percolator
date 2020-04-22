# - Try to find XSD
# Once done this will define
#
# This module defines
# XSD_INCLUDE_DIR, where to find elements.hxx, etc.
# XSD_EXECUTABLE, the exe file
# XSD_FOUND, If false, don't try to use it.

FIND_PATH(XSD_INCLUDE_DIR xsd/cxx/parser/elements.hxx
PATHS
  $ENV{XSDDIR}
  /usr/local/opt
  /usr/local
  /usr
PATH_SUFFIXES
  include
  libxsd
  xsd
  NO_CMAKE_FIND_ROOT_PATH
)

FIND_PROGRAM(XSD_EXECUTABLE
  NAMES
    xsdcxx xsd
  PATHS
    $ENV{XSDDIR}/bin
    /usr/local/Cellar/bin
    /usr/local/bin
    /usr/bin
  NO_SYSTEM_ENVIRONMENT_PATH
)


# if the include and the program are found then we have it
IF(XSD_INCLUDE_DIR)
  IF(XSD_EXECUTABLE)
    SET( XSD_FOUND "YES" )
  ENDIF(XSD_EXECUTABLE)
ENDIF(XSD_INCLUDE_DIR)


MARK_AS_ADVANCED(
  XSD_INCLUDE_DIR
  XSD_EXECUTABLE
)
