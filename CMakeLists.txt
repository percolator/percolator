#############################################################################
# NAME AND VERSION
###############################################################################

cmake_minimum_required (VERSION 3.15)

if(NOT CRUX AND APPLE AND ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER 3.0)
  # Identify AppleClang's COMPILER_ID as AppleClang instead of just Clang
  # (https://gitlab.kitware.com/cmake/cmake/issues/19067)
  cmake_policy(SET CMP0025 NEW)
endif()

# The macOS binaries run on macOS 12 and later, not only on the version of
# the build machine. This has to be set before project(). A value given with
# -DCMAKE_OSX_DEPLOYMENT_TARGET or in the MACOSX_DEPLOYMENT_TARGET
# environment variable is kept.
if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin"
   AND NOT DEFINED CMAKE_OSX_DEPLOYMENT_TARGET
   AND "$ENV{MACOSX_DEPLOYMENT_TARGET}" STREQUAL "")
  set(CMAKE_OSX_DEPLOYMENT_TARGET "12.0" CACHE STRING
      "Minimum macOS version of the binaries")
endif()

project(PERCOLATOR)

# The oldest supported GCC is 9: the code uses std::filesystem, which GCC 8
# keeps in a separate library (the link fails with undefined references).
# GCC 11 adds std::to_chars and std::from_chars for double; with GCC 9 and 10
# the slower number conversions are used (see src/CMakeLists.txt).
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
  message(FATAL_ERROR "Percolator needs GCC 9 or later "
                      "(this is GCC ${CMAKE_CXX_COMPILER_VERSION}).")
endif()
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
get_arch()

# STORE NEWLY SET VARIABLES IN *.h.cmake FILES
load_config_files("${CMAKE_CURRENT_SOURCE_DIR}")

###############################################################################
# COMPILING
###############################################################################

# MSVC: link the C and C++ runtime statically (/MT), so that percolator.exe
# and qvality.exe run without the Visual C++ Redistributable. This has to be
# set before the targets (including GoogleTest's) are defined. OpenMP's
# vcomp140.dll has no static version; it is installed next to the
# executables (see INSTALLING below).
option(USE_MSVC_STATIC_RUNTIME "Link the MSVC runtime statically (/MT)" ON)
if(MSVC AND USE_MSVC_STATIC_RUNTIME AND NOT DEFINED CMAKE_MSVC_RUNTIME_LIBRARY)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

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


# MSVC: install the redistributable DLLs, among them vcomp140.dll (OpenMP),
# next to the executables, which makes the installed bin directory portable.
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

include(CPack)
