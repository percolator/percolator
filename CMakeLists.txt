###############################################################################
# PREPARING TO INSTALL
###############################################################################
macro(my_set variable default options)
  if(NOT DEFINED ${variable})
    set(${variable} ${default} CACHE STRING ${options} FORCE)
  endif(NOT DEFINED ${variable})
endmacro(my_set)

# CMAKE_INSTALL_PREFIX must be set from command line with -D 
my_set(GOOGLE_TEST "FALSE" "Select if Google test should be included, options are: TRUE or FALSE") 
my_set(GOOGLE_TEST_PATH "../gtest-1.5.0" "Give path to google test code")
my_set(CMAKE_BUILD_TYPE "Release" "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
my_set(CMAKE_PREFIX_PATH "../" "Default path to packages")

# PRINT VARIBALES TO STDOUT
MESSAGE( STATUS )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS " Building PERCOLATOR:" )
MESSAGE( STATUS "change a configuration variable with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}" )
MESSAGE( STATUS "GOOGLE_TEST = ${GOOGLE_TEST}" )
MESSAGE( STATUS "GOOGLE_TEST_PATH = ${GOOGLE_TEST_PATH}" )
MESSAGE( STATUS "TARGET_ARCH = ${TARGET_ARCH}" )
MESSAGE( STATUS "TOOL CHAIN FILE = ${CMAKE_TOOLCHAIN_FILE}")
MESSAGE( STATUS "PROFILING = ${PROFILING}")
#MESSAGE( STATUS "LOG_FILE = ${LOG_ERROR_FILE}")
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS )


###############################################################################
# NAME AND VERSION
###############################################################################

cmake_minimum_required (VERSION 2.6)
project(PERCOLATOR)
 
# VERSION AND NAMESPACES
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "00")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(PIN_VERSION_MAJOR "1")
set(PIN_VERSION_MINOR "3")
set(POUT_VERSION_MAJOR "1")
set(POUT_VERSION_MINOR "4")
set(LOG_FILE ${LOG_ERROR_FILE})

IF(WIN32) 
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".dll" ${CMAKE_FIND_LIBRARY_SUFFIXES})
    SET(CMAKE_FIND_LIBRARY_PREFIXES "" "lib") 
    SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/lib/include ${CMAKE_SOURCE_DIR}/lib/include/log4cxx)
    SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${CMAKE_SOURCE_DIR}/lib/lib)
ELSEIF(UNIX)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so" ".la" ${CMAKE_FIND_LIBRARY_SUFFIXES})
    SET(CMAKE_FIND_LIBRARY_PREFIXES "" "lib")
    SET(CMAKE_INCLUDE_PATH /usr/include ${CMAKE_INCLUDE_PATH})
    SET(CMAKE_LIBRARY_PATH /usr/lib /usr/lib64 ${CMAKE_LIBRARY_PATH})
ENDIF(WIN32)

# Detect machine architecture, on UNIX:
if(UNIX AND NOT APPLE AND NOT DEFINED ${TARGET_ARCH})
	execute_process(COMMAND dpkg --print-architecture
		OUTPUT_VARIABLE TARGET_ARCH
                OUTPUT_STRIP_TRAILING_WHITESPACE)
	message(STATUS "Architecture: ${TARGET_ARCH}")
endif(UNIX AND NOT APPLE AND NOT DEFINED ${TARGET_ARCH})
if(APPLE AND NOT DEFINED ${TARGET_ARCH})
	execute_process(COMMAND uname -m
		OUTPUT_VARIABLE TARGET_ARCH
                OUTPUT_STRIP_TRAILING_WHITESPACE)
	message(STATUS "Architecture: ${TARGET_ARCH}")
endif(APPLE AND NOT DEFINED ${TARGET_ARCH})

if(MINGW OR WIN32) # windows
  set(WRITABLE_DIR "C:\\\\Docume~1\\\\AllUse~1\\\\Documents\\\\percolator\\\\")
  set(PIN_SCHEMA_LOCATION "${WRITABLE_DIR}xml-pin-${PIN_VERSION_MAJOR}-${PIN_VERSION_MINOR}\\\\")
  set(POUT_SCHEMA_LOCATION "${WRITABLE_DIR}xml-pout-${POUT_VERSION_MAJOR}-${POUT_VERSION_MINOR}\\\\")
  set(TEMP_DIR "%TEMP%\\\\")
  add_definitions(-D_WIN32)
else() # linux
  set(WRITABLE_DIR "${CMAKE_INSTALL_PREFIX}/share/percolator/")
  set(PIN_SCHEMA_LOCATION "${WRITABLE_DIR}xml-pin-${PIN_VERSION_MAJOR}-${PIN_VERSION_MINOR}/")
  set(POUT_SCHEMA_LOCATION "${WRITABLE_DIR}xml-pout-${POUT_VERSION_MAJOR}-${POUT_VERSION_MINOR}/")
  set(TEMP_DIR "/tmp/")
endif()
set(PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(percolator-in-namespace http://per-colator.com/percolator_in/${PIN_VERSION_MAJOR}${PIN_VERSION_MINOR})
set(percolator-out-namespace http://per-colator.com/percolator_out/${POUT_VERSION_MAJOR}${POUT_VERSION_MINOR})

# STORE NEWLY SET VARIABLES IN *.h.cmake FILES
file(GLOB_RECURSE configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/${myReplFile} @ONLY)
endforeach()


###############################################################################
# COMPILING
###############################################################################

# STATIC/DYNAMIC LIBRARY INCLUSION
if(MINGW OR WIN32 OR APPLE)
  OPTION(STATIC "build a static binaries. This is mainly intended for the win32-package and the MACOSX-bundle." ON )
else(MINGW OR WIN32 OR APPLE)
  OPTION(STATIC "build a static binaries. This is mainly intended for the win32-package and the MACOSX-bundle." OFF )
endif()

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

###############################################################################
# TESTING
###############################################################################

# Enabling system level tests (Ctest)
enable_testing()
# Scheduling system level tests
add_subdirectory(data/system_tests/percolator)
# Scheduling unit level tests
if(GOOGLE_TEST)
  add_subdirectory(data/unit_tests/percolator)
endif()

###############################################################################
# INSTALLING
###############################################################################

if(NOT MINGW)
  install(FILES src/xml/percolator_in.xsd DESTINATION ${PIN_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
  install(FILES src/xml/percolator_out.xsd DESTINATION ${POUT_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
endif()

###############################################################################
# PACKAGING
###############################################################################

# PACKAGING OPTIONS: GENERAL
set(PACKAGE_NICKNAME "percolator")
set(CPACK_SET_DESTDIR TRUE) # Enabling absolute paths for CPack (important!)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Percolator and qvality, two post processors for shotgun proteomics data.")
set(CPACK_SOURCE_GENERATOR "TGZ") # This file format is used to package source code ("make package_source")
set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_CURRENT_BINARY_DIR};${CMAKE_PROJECT_NAME};ALL;/")
set(CPACK_PACKAGE_VENDOR "Stockholms Universitet")
set(CPACK_PACKAGE_NAME ${PACKAGE_NICKNAME})
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PERCOLATOR_SOURCE_DIR}/src/README")
set(CPACK_RESOURCE_FILE_LICENSE "${PERCOLATOR_SOURCE_DIR}/src/COPYING")
set(PERLUQVAL_VERSION "v${CPACK_PACKAGE_VERSION_MAJOR}_${CPACK_PACKAGE_VERSION_MINOR}")
my_set(TARGET_ARCH "noarch" "set build architecture: options are i686, i386 and amd64")

IF(APPLE)
	set(CPACK_GENERATOR "TGZ;TBZ2;OSXX11")
        set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-mac-${TARGET_ARCH}")
        set(CPACK_SOURCE_PACKAGE_FILE_NAME  "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-mac-src")
        if(${TARGET_ARCH} STREQUAL "x86_64")
	  set(CMAKE_OSX_ARCHITECTURES "ppc;i386;ppc64;x86_64" CACHE STRING "Build architectures for OSX" FORCE)
        else(${TARGET_ARCH} STREQUAL "x86_64")
	  set(CMAKE_OSX_ARCHITECTURES ppc;i386 CACHE STRING "Build architectures for OSX" FORCE)
        endif(${TARGET_ARCH} STREQUAL "x86_64")
        set(MACOSX_BUNDLE_BUNDLE_NAME Converters)
        set(MACOSX_BUNDLE_BUNDLE_VERSION ${VERSION})
        set(MACOSX_BUNDLE_SHORT_VERSION_STRING ${VERSION})
        set(MACOSX_BUNDLE_LONG_VERSION_STRING "Version ${VERSION}")
elseif(UNIX)
	set(CPACK_GENERATOR "DEB;RPM") # Distribution packaging formats ("make package"); RPM packaging can also be included later
        set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-linux-${TARGET_ARCH}")
        set(CPACK_SOURCE_PACKAGE_FILE_NAME  "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-linux-src")
endif()

set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}")
set(CPACK_STRIP_FILES TRUE)

file(GLOB_RECURSE DOT_FILES_BEGIN ".*") # To be ignored by Cpack
file(GLOB_RECURSE TILD_FILES "*~*") # To be ignored by Cpack

set(CPACK_SOURCE_IGNORE_FILES "/CVS/;/.svn/;/.swp$/;cscope.*;/.git/;${CMAKE_CURRENT_BINARY_DIR}/;/.bzr/;/.settings/;${DOT_FILES_BEGIN};${TILD_FILES}")

# PACKAGING OPTIONS: DEB
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Lukas Käll <lukas.kall@cbr.su.se>")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${TARGET_ARCH})
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxerces-c-dev, libboost-dev, libboost-filesystem-dev (>=1.46), libboost-system-dev, libboost-thread-dev")

# PACKAGING OPTIONS: RPM
set(CPACK_RPM_PACKAGE_LICENSE "Apache 2.0 license")
set(CPACK_RPM_PACKAGE_VENDOR "Stockholms Universitet")
set(CPACK_RPM_PACKAGE_GROUP "Käll-lab")
set(CPACK_RPM_PACKAGE_DEPENDS "libc6, libgcc1")
set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS}, xerces-c-devel, boost-devel, boost-filesystem (>=1.46), boost-system, boost-thread")
set(CPACK_RPM_PACKAGE_RELOCATABLE FALSE)

# NSIS package options
# Requires NSIS for creating the package 
if(MINGW)
  add_custom_command( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/percolator-${PACKAGE_VERSION}-win32.exe" COMMAND makensis "${CMAKE_BINARY_DIR}/cmake/windows_percolator/script.nsi" DEPENDS percolator qvality) 
  add_custom_target( win32installer DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/percolator-${PACKAGE_VERSION}-win32.exe" )
endif()

IF(MINGW)

  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin")
  file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
             FILES_MATCHING PATTERN "boost_system*.dll")
  file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
             FILES_MATCHING PATTERN "boost_filesystem*.dll")

  file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
             FILES_MATCHING PATTERN "libxerces-c*.dll")
             
  #file(COPY "${MING_PATH}/bin/"
             #DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
             #FILES_MATCHING PATTERN "boost_thread*.dll")

  file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
             FILES_MATCHING PATTERN "libgcc_s_sjlj*.dll")

  file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
             FILES_MATCHING PATTERN "libstdc++*.dll")

endif()

include(CPack)