#############################################################################
# NAME AND VERSION
###############################################################################
cmake_minimum_required (VERSION 2.6)

if(NOT CRUX AND APPLE AND ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER 3.0)
  # Identify AppleClang's COMPILER_ID as AppleClang instead of just Clang 
  # (https://gitlab.kitware.com/cmake/cmake/issues/19067)
  cmake_policy(SET CMP0025 NEW)
endif()

project(PERCOLATOR)
set(PERCOLATOR_SOURCE_DIR ${CMAKE_SOURCE_DIR})
include("CommonCMake.txt")

###############################################################################
# PREPARING TO INSTALL
###############################################################################

my_set(CMAKE_BUILD_TYPE "Debug" "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
my_set(CMAKE_PREFIX_PATH "../" "Default path to packages")
option(XML_SUPPORT "Choose to support xml input (slower compilation)." OFF)
if(XML_SUPPORT)
  add_definitions(-DXML_SUPPORT)
endif(XML_SUPPORT)

# PRINT VARIBALES TO STDOUT
MESSAGE( STATUS )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS "Building PERCOLATOR:" )
MESSAGE( STATUS "change a configuration variable with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}" )
MESSAGE( STATUS "XML_SUPPORT = ${XML_SUPPORT}" )
MESSAGE( STATUS "GOOGLE_TEST = ${GOOGLE_TEST}" )
MESSAGE( STATUS "GOOGLE_TEST_PATH = ${GOOGLE_TEST_PATH}" )
MESSAGE( STATUS "TARGET_ARCH = ${TARGET_ARCH}" )
MESSAGE( STATUS "TOOL CHAIN FILE = ${CMAKE_TOOLCHAIN_FILE}")
MESSAGE( STATUS "PROFILING = ${PROFILING}")
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS )


# Detect machine architecture, on UNIX:
get_arch("${TARGET_ARCH}")

# STORE NEWLY SET VARIABLES IN *.h.cmake FILES
load_config_files("${CMAKE_CURRENT_SOURCE_DIR}")

###############################################################################
# COMPILING
###############################################################################
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


if( MINGW AND XML_SUPPORT )
    message( STATUS "    Installing system-libraries: MinGW DLLs." )
#    set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${MINGW_PATH}/bin/libgcc_s_sjlj-1.dll ${MINGW_PATH}/bin/libstdc++-6.dll ${MINGW_PATH}/bin/libxerces-c-3-1.dll ${MINGW_PATH}/bin/libboost_system.dll ${MINGW_PATH}/bin/libboost_filesystem.dll)
    set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${MINGW_PATH}/bin/libxerces-c-3-1.dll)
    include( InstallRequiredSystemLibraries )
endif( MINGW AND XML_SUPPORT )

if( MSVC )
  set( CMAKE_INSTALL_OPENMP_LIBRARIES TRUE )
  if ( XML_SUPPORT )
    message( STATUS "    Installing system-libraries: Win32 DLLs." )
    # TODO: Some hacking is needed to install the .dll instead of the .lib. CMake plans to include these DLL variables in the future
    string(REPLACE "xerces-c_3.lib" "../bin/xerces-c_3_1.dll;" XERCESC_DLL_LIBRARIES ${XERCESC_LIBRARIES})
    message( STATUS "      Installing Xerces DLL: ${XERCESC_DLL_LIBRARIES}")
    set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${XERCESC_DLL_LIBRARIES} )
  endif ( XML_SUPPORT )
  include( InstallRequiredSystemLibraries )
endif( MSVC )

install(FILES src/xml/percolator_in.xsd DESTINATION ${XML_REL_PATH}/${PIN_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
install(FILES src/xml/percolator_out.xsd DESTINATION ${XML_REL_PATH}/${POUT_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)

# ###############################################################################
# # PACKAGING
# ###############################################################################
# # PACKAGING OPTIONS: GENERAL
if(NOT XML_SUPPORT)
  set(NICKNAME_SUFFIX "-noxml")
endif()
set(PACKAGE_NICKNAME "percolator${NICKNAME_SUFFIX}")

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Percolator and qvality, two post processors for shotgun proteomics data.")
include("CPack.txt")

set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS} libc6, libgcc1")



include(CPack)
