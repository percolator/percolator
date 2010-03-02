# The name of the project is "LUMINATOR". CMakeLists files in this project can
# refer to the root source directory of the project as ${LUMINATOR_SOURCE_DIR} and
# to the root binary directory of the project as ${LUMINATOR_BINARY_DIR}.
cmake_minimum_required(VERSION 2.6)
project(LUMINATOR)

# Specifiying directory where built executables will be placed
set(EXECUTABLE_OUTPUT_PATH "${LUMINATOR_SOURCE_DIR}/bin")

# Going down to sources directory
add_subdirectory(src)

# Packaging general options
set(PACKAGE_NICKNAME "perluqval")
set(CPACK_SET_DESTDIR TRUE) # Enabling absolute paths for CPack (important!)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Luminator and Qvalumity software package for protein classification using unique algorithms.")
set(CPACK_GENERATOR "TGZ;DEB") # Distribution packaging formats ("make package"); RPM packaging can also be included later
set(CPACK_SOURCE_GENERATOR "TGZ") # This file format is used to package source code ("make package_source")
set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_CURRENT_BINARY_DIR};${CMAKE_PROJECT_NAME};ALL;/")
set(CPACK_PACKAGE_VENDOR "Matrix Science Ltd.")
set(CPACK_PACKAGE_NAME ${PACKAGE_NICKNAME})
set(CPACK_PACKAGE_DESCRIPTION_FILE "${LUMINATOR_SOURCE_DIR}/src/README")
set(CPACK_RESOURCE_FILE_LICENSE "${LUMINATOR_SOURCE_DIR}/src/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "1x")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(PERLUQVAL_VERSION "v${CPACK_PACKAGE_VERSION_MAJOR}_${CPACK_PACKAGE_VERSION_MINOR}")
set(TARGET_ARCH "amd64") # Suitable platform architecture for distribution package
set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-linux-${TARGET_ARCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME  "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-src")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}")
set(CPACK_STRIP_FILES TRUE)
set(CPACK_SOURCE_IGNORE_FILES "/CVS/;/.svn/;/.swp$/;.*~;cscope.*;/.git/;${CMAKE_CURRENT_BINARY_DIR}/;/.bzr/;/.settings/;.cdt*;.proj*;.cproj*;.cvs*;~$;.gitig*")

# DEB package options
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Lukas Käll <lukas.kall@cbr.su.se>, Luminita Moruz <lumi.moruz@gmail.com>")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${TARGET_ARCH})
#set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6, libgcc1") # Actual dependencies should be discovered later (TODO)

# RPM package options
set(CPACK_RPM_PACKAGE_LICENSE "Apache 2.0 license")
set(CPACK_RPM_PACKAGE_VENDOR "Matrix Science Ltd.")
set(CPACK_RPM_PACKAGE_GROUP "Lukas Käll group")
#set(CPACK_RPM_PACKAGE_REQUIRES "libc6, libgcc1") # Actual dependencies should be discovered later (TODO)

# Use CPack module within CMake
include(CPack)
