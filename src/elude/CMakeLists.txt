cmake_minimum_required(VERSION 2.8)
project(ELUDE)

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

#TODO: VERSIONING??
# VERSION AND NAMESPACES
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "00")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

if(STATIC AND MINGW) # windows
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
if(STATIC AND MINGW) # windows
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

target_link_libraries(elude eludelibrary)

#TODO: HERE THERE SHOULD BE STUFF ABOUT LIBGCC

IF(STATIC)
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
# TESTING
###############################################################################

#TODO: make sure this works with cpack
if(NOT MINGW)
  install(TARGETS elude EXPORT ELUDE DESTINATION ./bin) # Important to use relative path here (used by CPack)!
  install(DIRECTORY data/elude/models DESTINATION ${ELUDE_MODELS_PATH}.. FILE_PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
endif()

###############################################################################
# PACKAGING
###############################################################################

# PACKAGING OPTIONS: GENERAL
set(PACKAGE_NICKNAME "elude")
set(CPACK_SET_DESTDIR TRUE) # Enabling absolute paths for CPack (important!)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Elude software package for prediction of retention times in mass spectorometry experiments.")
set(CPACK_GENERATOR "DEB;RPM") # Distribution packaging formats ("make package"); RPM packaging can also be included later
set(CPACK_SOURCE_GENERATOR "TGZ") # This file format is used to package source code ("make package_source")
set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_CURRENT_BINARY_DIR};${CMAKE_PROJECT_NAME};ALL;/")
set(CPACK_PACKAGE_VENDOR "Stockholms Universitet")
set(CPACK_PACKAGE_NAME ${PACKAGE_NICKNAME})
#TODO: ADD README FILE
#set(CPACK_PACKAGE_DESCRIPTION_FILE "${PERCOLATOR_SOURCE_DIR}/src/elude/README")
set(CPACK_RESOURCE_FILE_LICENSE "${PERCOLATOR_SOURCE_DIR}/src/elude/COPYING")
set(PERLUQVAL_VERSION "v${CPACK_PACKAGE_VERSION_MAJOR}_${CPACK_PACKAGE_VERSION_MINOR}")
my_set(TARGET_ARCH "amd64" "set build architecture: options are i386 and amd64")
set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-linux-${TARGET_ARCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME  "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-src")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}")
set(CPACK_STRIP_FILES TRUE)#       NOT SURE IF THIS IS RECURSIVE...

file(GLOB_RECURSE DOT_FILES_BEGIN ".*") # To be ignored by Cpack
file(GLOB_RECURSE TILD_FILES "*~*") # To be ignored by Cpack
# uncomment to print list of files ignored by Cpack
#message("ignoring the following:")
#message("${DOT_FILES_BEGIN}")
#message("${TILD_FILES}")
set(CPACK_SOURCE_IGNORE_FILES "/CVS/;/.svn/;/.swp$/;cscope.*;/.git/;${CMAKE_CURRENT_BINARY_DIR}/;/.bzr/;/.settings/;${DOT_FILES_BEGIN};${TILD_FILES}")

#TODO: I DONT THINK ELUDE NEEDS ANY LLIBRARIES, BUT DOUBLE CHECK...
# PACKAGING OPTIONS: DEB
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Lukas KÃ€ll <lukas.kall@cbr.su.se>, Luminita Moruz <lumi.moruz@gmail.com>")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${TARGET_ARCH})
#set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxerces-c-dev, libboost-dev")

#TODO: I DONT THINK ELUDE NEEDS ANY LLIBRARIES, BUT DOUBLE CHECK...
# PACKAGING OPTIONS: RPM
set(CPACK_RPM_PACKAGE_LICENSE "Apache 2.0 license")
set(CPACK_RPM_PACKAGE_VENDOR "Stockholms Universitet")
set(CPACK_RPM_PACKAGE_GROUP "Kll-lab")
set(CPACK_RPM_PACKAGE_DEPENDS "libc6, libgcc1")
#set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS}, xerces-c-devel, boost-devel")

# NSIS package options
# Requires NSIS for creating the package 
add_custom_command( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/elude-${PACKAGE_VERSION}-win32.exe" COMMAND makensis "${CMAKE_BINARY_DIR}/cmake/windows_elude/script.nsi" DEPENDS elude) 
add_custom_target( win32installer DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/elude-${PACKAGE_VERSION}-win32.exe" )

# Use CPack
# the following MUST appear AFTER CPACK variables have been set
include(CPack)
