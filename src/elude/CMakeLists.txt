cmake_minimum_required (VERSION 2.8)
project(PERCOLATOR_ELUDE)

set(PERCOLATOR_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")

# PRINT VARIBALES
MESSAGE( STATUS )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS "Building ELUDE:" )
MESSAGE( STATUS "change a configuration variable with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "PERCOLATOR_SOURCE_DIR = ${PERCOLATOR_SOURCE_DIR}" )
MESSAGE( STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}" )
MESSAGE( STATUS "GOOGLE_TEST_PATH = ${GOOGLE_TEST_PATH}" )
MESSAGE( STATUS "TARGET_ARCH = ${TARGET_ARCH}" )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS )

# VERSION AND NAMESPACES
set(CPACK_PACKAGE_VERSION_MAJOR "2")
set(CPACK_PACKAGE_VERSION_MINOR "02")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(PIN_VERSION_MAJOR "1")
set(PIN_VERSION_MINOR "2")
set(POUT_VERSION_MAJOR "1")
set(POUT_VERSION_MINOR "3")
set(PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

# Detect machine architecture, on UNIX:
if(UNIX)
	execute_process(COMMAND uname -m
		OUTPUT_VARIABLE TARGET_ARCH
                OUTPUT_STRIP_TRAILING_WHITESPACE)
	message(STATUS "Architecture uname -m: ${TARGET_ARCH}")
endif(UNIX)


if(MINGW OR WIN32) # windows
  set(ELUDE_WRITABLE_DIR "C:\\\\Docume~1\\\\AllUse~1\\\\Documents\\\\elude\\\\")
  set(ELUDE_MODELS_PATH "C:\\\\Docume~1\\\\AllUse~1\\\\Documents\\\\elude\\\\models\\\\")
  set(TEMP_DIR "%TEMP%\\\\")
  add_definitions(-D_WIN32)
else() # linux
  set(ELUDE_WRITABLE_DIR "${CMAKE_INSTALL_PREFIX}/share/elude/")
  set(ELUDE_MODELS_PATH ${CMAKE_INSTALL_PREFIX}/share/elude/models/)
  set(TEMP_DIR "/tmp/")
endif()

#SET USEFUL VARIABLES TO BE USED IN THE CODE
if(MINGW OR WIN32) # windows
  set(WRITABLE_DIR "C:\\\\Docume~1\\\\AllUse~1\\\\Documents\\\\percolator\\\\")
  set(PIN_SCHEMA_LOCATION "${WRITABLE_DIR}xml-pin-${PIN_VERSION_MAJOR}-${PIN_VERSION_MINOR}\\\\")
else() # linux
  set(WRITABLE_DIR "${CMAKE_INSTALL_PREFIX}/share/percolator/")
  set(PIN_SCHEMA_LOCATION "${WRITABLE_DIR}xml-pin-${PIN_VERSION_MAJOR}-${PIN_VERSION_MINOR}/")
endif()

set(PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(percolator-in-namespace http://per-colator.com/percolator_in/${PIN_VERSION_MAJOR}${PIN_VERSION_MINOR})


# STORE NEWLY SET VARIABLES IN *.h.cmake FILES
file(GLOB configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/../ ${CMAKE_CURRENT_SOURCE_DIR}/../*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/${myReplFile} @ONLY)
endforeach()

file(GLOB configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/windows_elude ${CMAKE_CURRENT_SOURCE_DIR}/cmake/windows_elude/*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/windows_elude/${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/cmake/windows_elude/${myReplFile} @ONLY)
endforeach()


## FIND BOOST PACKAGE
SET (Boost_ADDITIONAL_VERSIONS "1.36" "1.36.0" "1.41" "1.41.0" "1.39" "1.39.0" "1.42.0" "1.42" "1.43.0" 
			       "1.43." "1.44.0" "1.44.0" "1.45.0" "1.45" "1.46.0" "1.46" "1.47.0" "1.47")
find_package(Boost REQUIRED)
if(Boost_FOUND)
  message(STATUS "Package Boost-${Boost_VERSION} found : ${Boost_INCLUDE_DIR}")
else(Boost_FOUND)
  message(FATAL_ERROR "The package Boost has not been found")
endif(Boost_FOUND)  
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_INCLUDE_DIRS})


###############################################################################
# COMPILING
###############################################################################
# include and link code from percolator
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${PERCOLATOR_SOURCE_DIR}/src)
link_directories(${PERCOLATOR_SOURCE_DIR}/src)

file(GLOB ELUDE_SOURCES *.cpp)
file(GLOB P_SOURCES ${PERCOLATOR_SOURCE_DIR}/src/PSMDescription.cpp)

add_library(eludelibrary RetentionFeatures.cpp DataManager.cpp EludeMain.cpp LibSVRModel.cpp LibsvmWrapper.cpp SVRModel.h RetentionModel.cpp EludeCaller.cpp  LTSRegression.cpp ../svm.cpp ../Normalizer.cpp ../UniNormalizer.cpp ../StdvNormalizer.cpp ../Option.cpp ../Enzyme.cpp ../PSMDescription.cpp ../Globals.cpp)

add_executable(elude ${ELUDE_SOURCES})

target_link_libraries(elude eludelibrary ${Boost_LIBRARIES})

#TODO: HERE THERE SHOULD BE STUFF ABOUT LIBGCC

IF(MINGW OR WIN32)
  set_property(TARGET elude PROPERTY LINK_SEARCH_END_STATIC ON)
  set_target_properties(elude PROPERTIES LINK_FLAGS "-Wl,-Bstatic -static-libgcc")
ENDIF()

###############################################################################
# TESTING
###############################################################################

# Enabling unit level tests (Ctest)
enable_testing()
# Scheduling unit level tests
if(GOOGLE_TEST)
  add_subdirectory(data/unit_tests/)
endif()
# Scheduling system level tests
add_subdirectory(data/system_tests/)


###############################################################################
# CROSS COMPILING
###############################################################################

#TODO: make sure the models are included in the windows installer
if(NOT MINGW)
  install(TARGETS elude EXPORT ELUDE DESTINATION ./bin) # Important to use relative path here (used by CPack)!
  install(DIRECTORY data/elude/models DESTINATION ${ELUDE_MODELS_PATH}.. FILE_PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
endif()

###############################################################################
# PACKAGING
###############################################################################
macro(my_set variable default options)
  if(NOT DEFINED ${variable})
    set(${variable} ${default} CACHE STRING ${options} FORCE)
  endif(NOT DEFINED ${variable})
endmacro(my_set)

# PACKAGING OPTIONS: GENERAL
set(PACKAGE_NICKNAME "elude")
set(CPACK_SET_DESTDIR TRUE) # Enabling absolute paths for CPack (important!)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Elude software package for prediction of retention times in mass spectorometry experiments.")
set(CPACK_SOURCE_GENERATOR "TGZ") # This file format is used to package source code ("make package_source")
set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_CURRENT_BINARY_DIR};${CMAKE_PROJECT_NAME};ALL;/")
set(CPACK_PACKAGE_VENDOR "Stockholms Universitet")
set(CPACK_PACKAGE_NAME ${PACKAGE_NICKNAME})
set(CPACK_RESOURCE_FILE_LICENSE "${PERCOLATOR_SOURCE_DIR}/src/elude/COPYING")
set(PERLUQVAL_VERSION "v${CPACK_PACKAGE_VERSION_MAJOR}_${CPACK_PACKAGE_VERSION_MINOR}")
my_set(TARGET_ARCH "amd64" "set build architecture: options are i386 and amd64")
#TODO: ADD README FILE
IF(APPLE)
	SET(CPACK_GENERATOR "TGZ;TBZ2;OSXX11")
        set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-mac-${TARGET_ARCH}")
        set(CPACK_SOURCE_PACKAGE_FILE_NAME  "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-mac-src")
#        SET(CMAKE_OSX_SYSROOT /Developer/SDKs/MacOSX10.4u.sdk)
        SET(CMAKE_OSX_ARCHITECTURES ppc;i386)
#        SET(CMAKE_FIND_ROOT_PATH /Developer/SDKs/MacOSX10.4u.sdk)
        SET(MACOSX_BUNDLE_BUNDLE_NAME Elude)
        SET(MACOSX_BUNDLE_BUNDLE_VERSION ${VERSION})
        SET(MACOSX_BUNDLE_SHORT_VERSION_STRING ${VERSION})
        SET(MACOSX_BUNDLE_LONG_VERSION_STRING "Version ${VERSION}")
elseif(UNIX)
	set(CPACK_GENERATOR "DEB;RPM") # Distribution packaging formats ("make package"); RPM packaging can also be included later
        set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-linux-${TARGET_ARCH}")
        set(CPACK_SOURCE_PACKAGE_FILE_NAME  "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-mac-src")
endif()

set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}")
set(CPACK_STRIP_FILES TRUE)#       NOT SURE IF THIS IS RECURSIVE...

file(GLOB_RECURSE DOT_FILES_BEGIN ".*") # To be ignored by Cpack
file(GLOB_RECURSE TILD_FILES "*~*") # To be ignored by Cpack

set(CPACK_SOURCE_IGNORE_FILES "/CVS/;/.svn/;/.swp$/;cscope.*;/.git/;${CMAKE_CURRENT_BINARY_DIR}/;/.bzr/;/.settings/;${DOT_FILES_BEGIN};${TILD_FILES}")

# PACKAGING OPTIONS: DEB
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Lukas KÃ€ll <lukas.kall@cbr.su.se>, Luminita Moruz <lumi.moruz@gmail.com>")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${TARGET_ARCH})
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libboost-dev")

#TODO: I DONT THINK ELUDE NEEDS ANY LLIBRARIES, BUT DOUBLE CHECK...
# PACKAGING OPTIONS: RPM
set(CPACK_RPM_PACKAGE_LICENSE "Apache 2.0 license")
set(CPACK_RPM_PACKAGE_VENDOR "Stockholms Universitet")
set(CPACK_RPM_PACKAGE_GROUP "Kll-lab")
set(CPACK_RPM_PACKAGE_DEPENDS "libc6, libgcc1, boost-devel")

# NSIS package options
# Requires NSIS for creating the package 
if(MINGW)
  add_custom_command( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/elude-${PACKAGE_VERSION}-win32.exe" COMMAND makensis "${CMAKE_BINARY_DIR}/cmake/windows_elude/script.nsi" DEPENDS elude) 
  add_custom_target( win32installer DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/elude-${PACKAGE_VERSION}-win32.exe" )
endif()
# Use CPack
# the following MUST appear AFTER CPACK variables have been set
include(CPack)
