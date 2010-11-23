###############################################################################
# PREPARING TO INSTALL
###############################################################################

set(base "/tmp")
set(percolatorBuild "${base}/percolatorBuild")
# CMAKE_INSTALL_PREFIX must be set from command line with -D
# set(CMAKE_INSTALL_PREFIX "${base}") 
set(CMAKE_BUILD_TYPE "Debug")
set(GOOGLE_TEST "TRUE") 
set(EXCLUDE_CONVERTERS "FALSE")
set(EXCLUDE_ELUDE "FALSE")
set(CMAKE_PREFIX_PATH "${base}/xsd-3.3.0-x86_64-linux-gnu/")
set(GOOGLE_TEST_PATH "${base}/gtest-1.5.0")


# the Ubuntu package xsdcxx is as of 2010-10-21 version 3.2 and we want version 3.3
# so we download it from the Codesynthesis home page
execute_process(COMMAND file xsd-3.3.0-x86_64-linux-gnu WORKING_DIRECTORY ${base} RESULT_VARIABLE notExists OUTPUT_QUIET ERROR_QUIET)
if(notExists)
  execute_process(COMMAND wget http://codesynthesis.com/download/xsd/3.3/linux-gnu/x86_64/xsd-3.3.0-x86_64-linux-gnu.tar.bz2 WORKING_DIRECTORY ${base})
  execute_process(COMMAND tar xfj xsd-3.3.0-x86_64-linux-gnu.tar.bz2 WORKING_DIRECTORY ${base})
  execute_process(COMMAND rm xsd-3.3.0-x86_64-linux-gnu.tar.bz2 WORKING_DIRECTORY ${base})
endif()
# the Ubuntu package xsdcxx is as of 2010-10-21 version 1.3 and we want version 1.5
# so we download it from the Googletest home page
execute_process(COMMAND file gtest-1.5.0 WORKING_DIRECTORY ${base} RESULT_VARIABLE notExists OUTPUT_QUIET ERROR_QUIET)
if(notExists)
  execute_process(COMMAND wget http://googletest.googlecode.com/files/gtest-1.5.0.tar.gz WORKING_DIRECTORY ${base})
  execute_process(COMMAND tar xzf gtest-1.5.0.tar.gz WORKING_DIRECTORY ${base})
  execute_process(COMMAND rm gtest-1.5.0.tar.gz WORKING_DIRECTORY ${base})
  execute_process(COMMAND mkdir build WORKING_DIRECTORY ${base}/gtest-1.5.0)
  execute_process(COMMAND cmake ../ WORKING_DIRECTORY ${base}/gtest-1.5.0/build)
  execute_process(COMMAND make WORKING_DIRECTORY ${base}/gtest-1.5.0/build)
endif()

# PRINT VARIBALES TO STDOUT
MESSAGE( STATUS )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS "Building PERCOLATOR:" )
MESSAGE( STATUS "change a configuration variable with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "EXCLUDE_CONVERTERS = ${EXCLUDE_CONVERTERS}" )
MESSAGE( STATUS "EXCLUDE_ELUDE = ${EXCLUDE_ELUDE}" )
MESSAGE( STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}" )
MESSAGE( STATUS "GOOGLE_TEST_PATH = ${GOOGLE_TEST_PATH}" )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS )


###############################################################################
# NAME AND VERSION
###############################################################################

cmake_minimum_required(VERSION 2.8)
project(PERCOLATOR)

# VERSION AND NAMESPACES
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "17")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(SCHEMA_VERSION_MAJOR "1")
set(SCHEMA_VERSION_MINOR "1")

if(STATIC AND MINGW) # windows
  set(SCHEMA "C:\\\\Docume~1\\\\AllUse~1\\\\Documents\\\\percolator\\\\xml-${SCHEMA_VERSION_MAJOR}-${SCHEMA_VERSION_MINOR}\\\\")
  set(ELUDE_MODELS_PATH "C:\\\\Docume~1\\\\AllUse~1\\\\Documents\\\\percolator\\\\elude_model\\\\")
else() # linux
  set(ELUDE_MODELS_PATH ${CMAKE_INSTALL_PREFIX}/share/percolator/elude_model/)
  set(SCHEMA ${CMAKE_INSTALL_PREFIX}/share/percolator/xml-${SCHEMA_VERSION_MAJOR}-${SCHEMA_VERSION_MINOR}/)
endif()
set(PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(percolator-in-namespace http://per-colator.com/percolator_in/${SCHEMA_VERSION_MAJOR}${SCHEMA_VERSION_MINOR})
set(percolator-out-namespace http://per-colator.com/percolator_out/${SCHEMA_VERSION_MAJOR}${SCHEMA_VERSION_MINOR})

# STORE NEWLY SET VARIABLES IN *.h.cmake FILES
file(GLOB_RECURSE configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/${myReplFile} @ONLY)
endforeach()

###############################################################################
# COMPILING, INSTALLING AND TESTING
###############################################################################

# STATIC/DYNAMIC LIBRARY INCLUSION
OPTION(STATIC "build a static binaries. This is mainly intended for the win32-package and the MACOSX-bundle." OFF )

# useful macro
macro(my_find_program progname)
  find_program(${progname}_EXECUTABLE ${progname})
  if(NOT ${progname}_EXECUTABLE)
    message(FATAL_ERROR  "${progname} not found" )
  endif(NOT ${progname}_EXECUTABLE)
endmacro(my_find_program)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# COMPILE
add_subdirectory(src)

# Enabling unit level tests (Ctest)
enable_testing()
# Compiling unit level tests
add_subdirectory(data/tests)

#INSTALL
if(NOT MINGW)
  install(FILES src/xml/percolator_in.xsd DESTINATION ${SCHEMA} PERMISSIONS OWNER_READ WORLD_READ)
  install(FILES src/xml/percolator_out.xsd DESTINATION ${SCHEMA} PERMISSIONS OWNER_READ WORLD_READ)
  install(DIRECTORY data/elude_model DESTINATION ${SCHEMA}.. FILE_PERMISSIONS OWNER_READ OWNER_WRITE WORLD_READ WORLD_WRITE)
endif()

###############################################################################
# PACKAGING
###############################################################################

# PACKAGING OPTIONS: GENERAL
set(PACKAGE_NICKNAME "percolator")
set(CPACK_SET_DESTDIR TRUE) # Enabling absolute paths for CPack (important!)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Percolator, Elude and Qvality software package for protein classification using unique algorithms.")
set(CPACK_GENERATOR "DEB;RPM") # Distribution packaging formats ("make package"); RPM packaging can also be included later
set(CPACK_SOURCE_GENERATOR "TGZ") # This file format is used to package source code ("make package_source")
set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_CURRENT_BINARY_DIR};${CMAKE_PROJECT_NAME};ALL;/")
set(CPACK_PACKAGE_VENDOR "Stockholms Universitet")
set(CPACK_PACKAGE_NAME ${PACKAGE_NICKNAME})
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PERCOLATOR_SOURCE_DIR}/src/README")
set(CPACK_RESOURCE_FILE_LICENSE "${PERCOLATOR_SOURCE_DIR}/src/COPYING")
set(PERLUQVAL_VERSION "v${CPACK_PACKAGE_VERSION_MAJOR}_${CPACK_PACKAGE_VERSION_MINOR}")
set(TARGET_ARCH "amd64") # Suitable platform architecture for distribution package: 32bit=i386, 64bit=amd64
set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-linux-${TARGET_ARCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME  "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-src")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}")
set(CPACK_STRIP_FILES TRUE)

file(GLOB_RECURSE DOT_FILES_BEGIN ".*") # To be ignored by Cpack
file(GLOB_RECURSE TILD_FILES "*~*") # To be ignored by Cpack
# uncomment to print list of files ignored by Cpack
#message("ignoring the following:")
#message("${DOT_FILES_BEGIN}")
#message("${TILD_FILES}")
set(CPACK_SOURCE_IGNORE_FILES "/CVS/;/.svn/;/.swp$/;cscope.*;/.git/;${CMAKE_CURRENT_BINARY_DIR}/;/.bzr/;/.settings/;${DOT_FILES_BEGIN};${TILD_FILES}")

# PACKAGING OPTIONS: DEB
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Lukas KÃ¤ll <lukas.kall@cbr.su.se>, Luminita Moruz <lumi.moruz@gmail.com>")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${TARGET_ARCH})
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxerces-c-dev, libboost-dev")
if(NOT EXCLUDE_CONVERTERS)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, gengetopt, libtokyocabinet-dev, zlib1g-dev, libsqlite3-dev")
endif()

# PACKAGING OPTIONS: RPM
set(CPACK_RPM_PACKAGE_LICENSE "Apache 2.0 license")
set(CPACK_RPM_PACKAGE_VENDOR "Stockholms Universitet")
set(CPACK_RPM_PACKAGE_GROUP "KŠll-lab")
set(CPACK_RPM_PACKAGE_DEPENDS "libc6, libgcc1")
set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS}, xerces-c-devel, boost-devel")
if(NOT EXCLUDE_CONVERTERS)
  set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS}, gengetopt, tokyocabinet-devel, zlib-devel, sqlite3-devel")
endif()

# NSIS package options
# Requires NSIS for creating the package 
add_custom_command( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/percolator-${PACKAGE_VERSION}-win32.exe" COMMAND makensis "${CMAKE_BINARY_DIR}/cmake/script.nsi" DEPENDS percolator elude qvality ) 
add_custom_target( win32installer DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/percolator-${PACKAGE_VERSION}-win32.exe" )

# Use CPack
# the following MUST appear AFTER CPACK variables have been set
include(CPack)

