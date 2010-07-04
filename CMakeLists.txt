cmake_minimum_required(VERSION 2.8)
project(PERCOLATOR)

macro(my_find_program progname)
  find_program(${progname}_EXECUTABLE ${progname})
  if(NOT ${progname}_EXECUTABLE)
    message(FATAL_ERROR  "${progname} not found" )
  endif(NOT ${progname}_EXECUTABLE)
endmacro(my_find_program)

OPTION(STATIC "build a static binaries. This is mainly intended for the win32-package and the MACOSX-bundle." OFF )

set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "14b")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(percolator-in-namespace http://per-colator.com/percolator_in/10)

file(GLOB_RECURSE configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/${myReplFile} @ONLY)
endforeach()

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)


add_subdirectory(src)

# Packaging general options
set(PACKAGE_NICKNAME "percolator")
set(CPACK_SET_DESTDIR TRUE) # Enabling absolute paths for CPack (important!)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Percolator, Elude and Qvality software package for protein classification using unique algorithms.")
set(CPACK_GENERATOR "TGZ;DEB") # Distribution packaging formats ("make package"); RPM packaging can also be included later
set(CPACK_SOURCE_GENERATOR "TGZ") # This file format is used to package source code ("make package_source")
set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_CURRENT_BINARY_DIR};${CMAKE_PROJECT_NAME};ALL;/")
set(CPACK_PACKAGE_VENDOR "Matrix Science Ltd.")
set(CPACK_PACKAGE_NAME ${PACKAGE_NICKNAME})
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PERCOLATOR_SOURCE_DIR}/src/README")
set(CPACK_RESOURCE_FILE_LICENSE "${PERCOLATOR_SOURCE_DIR}/src/COPYING")
set(PERLUQVAL_VERSION "v${CPACK_PACKAGE_VERSION_MAJOR}_${CPACK_PACKAGE_VERSION_MINOR}")
set(TARGET_ARCH "i386") # Suitable platform architecture for distribution package
set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-linux-${TARGET_ARCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME  "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-src")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}")
set(CPACK_STRIP_FILES TRUE)

file(GLOB_RECURSE DOT_FILES_BEGIN ".*")
file(GLOB_RECURSE TILD_FILES "*~*")
message("${DOT_FILES_BEGIN}")
message("${TILD_FILES}")

set(CPACK_SOURCE_IGNORE_FILES "/CVS/;/.svn/;/.swp$/;cscope.*;/.git/;${CMAKE_CURRENT_BINARY_DIR}/;/.bzr/;/.settings/;${DOT_FILES_BEGIN};${TILD_FILES}")

# DEB package options
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Lukas KÃ¤ll <lukas.kall@cbr.su.se>, Luminita Moruz <lumi.moruz@gmail.com>")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${TARGET_ARCH})
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6, libgcc1, libgtest0") # Actual dependencies should be discovered later (TODO)

# RPM package options
set(CPACK_RPM_PACKAGE_LICENSE "Apache 2.0 license")
set(CPACK_RPM_PACKAGE_VENDOR "Stockholms Universitet")
set(CPACK_RPM_PACKAGE_GROUP "KŠll-lab")
#set(CPACK_RPM_PACKAGE_REQUIRES "libc6, libgcc1") # Actual dependencies should be discovered later (TODO)

# NSIS package options
# Requires NSIS for creating the package 

add_custom_command( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/percolator-${PACKAGE_VERSION}-win32.exe" COMMAND makensis "${CMAKE_BINARY_DIR}/script.nsi" DEPENDS percolator elude qvality ) 
add_custom_target( win32installer DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/percolator-${PACKAGE_VERSION}-win32.exe" )



# Adding data directories for installation.
install(DIRECTORY data DESTINATION bin)
install(DIRECTORY models DESTINATION bin)

# Use CPack module within CMake
include(CPack)

# I comment out the testing to avoid some compilation error /Erik Sjolund 2010-05-07

# Enabling CTest functionality
#enable_testing()

# Running system and unit level tests
#add_subdirectory(tests)

