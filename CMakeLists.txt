#############################################################################
# NAME AND VERSION
###############################################################################

cmake_minimum_required (VERSION 3.15)

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

my_set(CMAKE_BUILD_TYPE "Release" "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
my_set(CMAKE_PREFIX_PATH "../" "Default path to packages")

option(NO_TESTS "Build system and integration tests" OFF)

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
MESSAGE( STATUS "NO_TESTS = ${NO_TESTS}" )
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

if(NOT NO_TESTS)
  # Enabling system level tests (Ctest)
  enable_testing()
  # Scheduling system level tests
  add_subdirectory(data/system_tests/percolator)
  add_subdirectory(tests/integration_tests/percolator)
  # Scheduling unit level tests
  if(NOT WITHOUT_GTEST)
    add_subdirectory(tests/unit_tests/percolator)
  endif()
endif()

###############################################################################
# INSTALLING
###############################################################################


if( MSVC )
  set( CMAKE_INSTALL_OPENMP_LIBRARIES TRUE )
  include( InstallRequiredSystemLibraries )
endif( MSVC )


# ###############################################################################
# # PACKAGING
# ###############################################################################
# # PACKAGING OPTIONS: GENERAL
set(PACKAGE_NICKNAME "percolator${NICKNAME_SUFFIX}")

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Percolator and qvality, two post processors for shotgun proteomics data.")
include("CPack.txt")

set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS} libc6, libgcc1")



include(CPack)
