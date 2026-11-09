

###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY PERCOLATOR
###############################################################################

# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(${CMAKE_CURRENT_BINARY_DIR})

MESSAGE( STATUS "Checking packages:" )

set(COMMON_LIBRARIES "")
find_package(OpenMP)
# Help AppleClang find Homebrew's libomp if not auto-detected. Conda builds
# turn this off, so that they never link to a library outside the environment.
option(HOMEBREW_LIBOMP "On macOS, fall back to Homebrew's libomp if OpenMP is not found" ON)
if(APPLE AND NOT OpenMP_CXX_FOUND AND HOMEBREW_LIBOMP)
  # Common Homebrew prefix on Apple Silicon
  set(_HOMEBREW_OMP_ROOT "/opt/homebrew/opt/libomp")
  if (EXISTS "${_HOMEBREW_OMP_ROOT}/include" AND EXISTS "${_HOMEBREW_OMP_ROOT}/lib")
    # Hint include/lib to CMake's OpenMP finder
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp" CACHE STRING "" FORCE)
    set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "" FORCE)
    set(OpenMP_omp_LIBRARY "${_HOMEBREW_OMP_ROOT}/lib/libomp.dylib" CACHE FILEPATH "" FORCE)
    include_directories("${_HOMEBREW_OMP_ROOT}/include")
    # Retry discovery
    find_package(OpenMP)
  endif()
endif()
if (OpenMP_CXX_FOUND)
    list(APPEND COMMON_LIBRARIES OpenMP::OpenMP_CXX)
endif()



# Make find_package honor UPPER-CASE <PKG>_ROOT vars like BOOST_ROOT (CMake 3.30+)
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 OLD)
endif()
if(POLICY CMP0169)
  cmake_policy(SET CMP0169 NEW)
endif()
if(POLICY CMP0177)
  cmake_policy(SET CMP0177 NEW)
endif()



# --- Boost configuration (platform-aware) ---
# Windows needs static Boost to match the CRT. Linux and macOS use shared libs
# (macOS Homebrew ships shared-only; Linux distros vary but shared is always present).
if(MSVC)
  set(Boost_USE_STATIC_LIBS ON)
else()
  set(Boost_USE_STATIC_LIBS OFF)
endif()
set(Boost_USE_MULTITHREADED ON)

# Prefer Boost's classic FindBoost over boost-cmake configs
set(Boost_NO_BOOST_CMAKE ON CACHE BOOL "" FORCE)

# Require the components you actually use
if(APPLE)
  find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS filesystem REQUIRED)
else()
  # Boost.System is header-only since 1.69, and 1.89 removed its stub
  # library, so it is linked only where it exists.
  find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS filesystem OPTIONAL_COMPONENTS system REQUIRED)
endif()

add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

if(WIN32)
  add_definitions(-DBOOST_ALL_NO_LIB) # disable autolinking in boost
endif()

include_directories(${Boost_INCLUDE_DIRS})
list(APPEND COMMON_LIBRARIES Boost::filesystem)
if(NOT APPLE AND TARGET Boost::system)
  list(APPEND COMMON_LIBRARIES Boost::system)
endif()

if(MSVC)
  add_definitions (/MP /D _CRT_SECURE_NO_WARNINGS /D _CRT_NONSTDC_NO_DEPRECATE)
endif(MSVC)
option(USE_MSVC_STATIC_RUNTIME "Link MSVC runtime statically (/MT)" ON)
if (MSVC AND USE_MSVC_STATIC_RUNTIME)
  foreach(flag_var CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
                   CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE)
    if(${flag_var} MATCHES "/MD")
      string(REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
    endif()
  endforeach()
  set(Boost_USE_STATIC_RUNTIME ON)  # match Boost to the CRT choice
endif()

find_package(Threads REQUIRED)
list(APPEND COMMON_LIBRARIES Threads::Threads)


#########################################
# COMPILE BLAS
#########################################
add_subdirectory(dblas)
set (DBLAS_LIBRARIES ${DBLAS_LIBRARIES} dblas)

###############################################################################
# PULL Eigen
###############################################################################

include(FetchContent)

# Avoid Eigen's conflicting 'blas' target
set(EIGEN_BUILD_BLAS OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)

if(NOT OpenMP_CXX_FOUND)
  message(STATUS "OpenMP not found; building Eigen without it")
endif()

# Avoid adding Eigen's blas directory (explicitly ignore it)
set(EIGEN_DISABLE_BLAS ON CACHE BOOL "Disable Eigen's BLAS build" FORCE)

# EXCLUDE_FROM_ALL (CMake 3.28 or later) keeps Eigen's headers and CMake
# files out of "make install" and the packages; only its headers are used.
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
  EXCLUDE_FROM_ALL
)

FetchContent_MakeAvailable(eigen)
include_directories(${eigen_SOURCE_DIR})

### Patch for avoiding _finite

set(EIGEN_MATHFUNCTIONS_FILE "${eigen_SOURCE_DIR}/Eigen/src/Core/MathFunctions.h")

# Do the patch if not already patched
file(READ "${EIGEN_MATHFUNCTIONS_FILE}" MATHFUNCS_CONTENT)

string(REGEX REPLACE
  "#elif EIGEN_USE_STD_FPCLASSIFY[ \r\n]+[^\r\n]*using std::isfinite;[ \r\n]+[^\r\n]*return isfinite[^\r\n]*;"
  "#elif EIGEN_USE_STD_FPCLASSIFY\n    return !std::isnan(x) && !std::isinf(x); // PATCHED to avoid std::isfinite on MSVC"
  MATHFUNCS_CONTENT_PATCHED
  "${MATHFUNCS_CONTENT}"
)

file(WRITE "${EIGEN_MATHFUNCTIONS_FILE}" "${MATHFUNCS_CONTENT_PATCHED}")
### End patch

if (MSVC)
  add_compile_definitions(
    _HAS_STD_BYTE=0     # disables std::byte so 'byte' in Windows headers is unambiguous
    NOMINMAX            # avoid macro min/max conflicts
    _WIN32_WINNT=0x0601 # silence that SDK warning (Windows 7 or higher)
  )
endif()

#############################################################################
# COMPILE PERCOLATOR
###############################################################################

add_library(perclibrary STATIC BaseSpline.cpp MassHandler.cpp ResultHolder.cpp PSMDescription.cpp MonotoneRegressor.cpp IsotonicPEP.cpp
								  SetHandler.cpp StdvNormalizer.cpp svm.cpp Caller.cpp CompositionSorter.cpp CrossValidation.cpp Enzyme.cpp Globals.cpp NoNormalizer.cpp Normalizer.cpp
								  SanityCheck.cpp UniNormalizer.cpp DataSet.cpp TabFileValidator.cpp FeatureNames.cpp LogisticRegression.cpp Numerical.cpp Option.cpp PosteriorEstimator.cpp
								  ProteinProbEstimator.cpp ProteinFDRestimator.cpp ScoreHolder.cpp Scores.cpp PseudoRandom.cpp Set.cpp SqtSanityCheck.cpp ssl.cpp PackedVector.cpp
								  PackedMatrix.cpp Logger.cpp MyException.cpp ProteinScoreHolder.cpp PickedProteinInterface.cpp FeatureMemoryPool.cpp Timer.cpp TmpDir.cpp ValidateTabFile.cpp Vector.cpp)

if(OpenMP_CXX_FOUND)
  target_link_libraries(perclibrary PUBLIC ${COMMON_LIBRARIES} OpenMP::OpenMP_CXX)
else()
  target_link_libraries(perclibrary PUBLIC ${COMMON_LIBRARIES})
  target_compile_definitions(perclibrary PUBLIC EIGEN_DONT_PARALLELIZE)
endif()

###############################################################################
# COMPILE INTERNAL LIBRARIES
###############################################################################

#COMPILE PICKED PROTEIN
add_subdirectory(picked_protein)

#INCLUDE PERCOLATOR HEADERS
include_directories(.)

#INCLUDE PICKED PROTEIN HEADERS FOR PERCLIBRARY
include_directories(picked_protein)

################################################################################

# COMPILE PERCOLATOR
if(CRUX)
  ADD_DEFINITIONS(-DCRUX)
  include_directories(${CRUX} ${CRUX}/src ${EXT_BINARY_DIR})
ELSE(CRUX)
  add_executable(percolator main.cpp)

  if(PROFILING AND UNIX)
    set_target_properties(percolator PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -DNDEBUG ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g -DNDEBUG ${PROFILING_COMMON_FLAGS}")
  endif()

  # SET LIBRARIES FOR PERCOLATOR
  list(APPEND COMMON_LIBRARIES perclibrary ${DBLAS_LIBRARIES} picked_protein ${CURL_LIBRARIES} ${RPC_LIBRARIES})
  if(NOT MSVC)
    if(NOT APPLE)
      list(APPEND COMMON_LIBRARIES stdc++)
    endif(NOT APPLE)
  endif(NOT MSVC)

  if(COVERAGE)
    target_link_libraries(percolator PRIVATE -fprofile-arcs ${COMMON_LIBRARIES})
  else(COVERAGE)
    target_link_libraries(percolator PRIVATE ${COMMON_LIBRARIES})
  endif(COVERAGE)

  # INSTALL PERCOLATOR
  install(TARGETS percolator EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!
ENDIF(CRUX)

###############################################################################
# COMPILE QVALITY
###############################################################################

add_subdirectory(qvality)
