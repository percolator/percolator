FILE(GLOB my_local "$ENV{HOME}/src/xsd*")
FILE(GLOB my_tmp "/tmp/xsd*")

FIND_PATH(XSD_INCLUDE_DIR xsd/cxx/parser/elements.hxx PATH_SUFFIXES libxsd PATHS
  "[HKEY_CURRENT_USER\\software\\xsd\\include]"
  "[HKEY_CURRENT_USER]\\xsd\\include]"
  ${my_tmp}
  ${my_local}
  /usr/local/include/xsd
  /usr/include/xsd
)

FIND_PROGRAM(XSD_EXECUTABLE 
  NAMES 
    xsd xsdcxx
  PATHS
    "[HKEY_CURRENT_USER\\xsd\\bin]"
    ${my_local}/bin
    ${my_local}
    ${my_tmp}/bin
    /usr/local/bin
    /usr/bin
)

#MESSAGE("my_local=${my_local}")


# if the include and the program are found then we have it
IF(XSD_INCLUDE_DIR)
  IF(XSD_EXECUTABLE)
    SET( XSD_FOUND "YES" )
  ENDIF(XSD_EXECUTABLE)
ENDIF(XSD_INCLUDE_DIR)

#MESSAGE("XSD_INCLUDE_DIR=${XSD_INCLUDE_DIR}")
#MESSAGE("XSD_EXECUTABLE=${XSD_EXECUTABLE}")

MARK_AS_ADVANCED(
  XSD_INCLUDE_DIR
  XSD_EXECUTABLE
) 
