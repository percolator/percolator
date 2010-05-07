cmake_minimum_required(VERSION 2.8)
project(PERCOLATOR)

macro(my_find_program progname)
  find_program(${progname}_EXECUTABLE ${progname})
  if(NOT ${progname}_EXECUTABLE)
    message(FATAL_ERROR  "${progname} not found" )
  endif(NOT ${progname}_EXECUTABLE)
endmacro(my_find_program)

set(percolator-in-namespace http://per-colator.com/percolator_in/10)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

add_subdirectory(src)

# Packaging general options
set(CPACK_SET_DESTDIR TRUE) # Enabling absolute paths for CPack (important!)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Percolator and Qvality proteins classification project")
set(CPACK_GENERATOR "TGZ;DEB") # Distribution packaging formats ("make package"); RPM packaging can also be included later
set(CPACK_SOURCE_GENERATOR "TGZ") # This file format is used to package source code ("make package_source")
set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_CURRENT_BINARY_DIR};${CMAKE_PROJECT_NAME};ALL;/")
set(CPACK_PACKAGE_VENDOR "Matrix Science Ltd.")
set(CPACK_PACKAGE_NAME "perqval")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PERCOLATOR_SOURCE_DIR}/src/README")
set(CPACK_RESOURCE_FILE_LICENSE "${PERCOLATOR_SOURCE_DIR}/src/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "14")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(PERQVAL_VERSION "v${CPACK_PACKAGE_VERSION_MAJOR}_${CPACK_PACKAGE_VERSION_MINOR}")
set(TARGET_ARCH "amd64") # Suitable platform architecture for distribution package
set(CPACK_PACKAGE_FILE_NAME "perqval-${PERQVAL_VERSION}-linux-${TARGET_ARCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME  "perqval-${PERQVAL_VERSION}-src")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "perqval-${PERQVAL_VERSION}")
set(CPACK_STRIP_FILES TRUE)
set(CPACK_SOURCE_IGNORE_FILES "/CVS/;/.svn/;/.swp$/;.*~;cscope.*;/.git/;${CMAKE_CURRENT_BINARY_DIR}/;/.bzr/;/.settings/;.cdt*;.proj*;.cproj*;.cvs*;~$;.gitig*")

# DEB package options
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Lukas Käll <lukas.kall@cbr.su.se>")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${TARGET_ARCH})
#set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6, libgcc1") # Actual dependencies should be discovered later (TODO)

# RPM package options
set(CPACK_RPM_PACKAGE_LICENSE "Apache 2.0 license")
set(CPACK_RPM_PACKAGE_VENDOR "Matrix Science Ltd.")
set(CPACK_RPM_PACKAGE_GROUP "Lukas Käll group")
#set(CPACK_RPM_PACKAGE_REQUIRES "libc6, libgcc1") # Actual dependencies should be discovered later (TODO)

# Use CPack module within CMake
include(CPack)

# I comment out the testing to avoid some compilation error /Erik Sjolund 2010-05-07

# Enabling CTest functionality
#enable_testing()

# Running system and unit level tests
#add_subdirectory(tests)

