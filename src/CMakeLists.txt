

###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY PERCOLATOR
###############################################################################

# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(${CMAKE_CURRENT_BINARY_DIR})

MESSAGE( STATUS "Checking packages:" )

set(COMMON_LIBRARIES "")
find_package(OpenMP)
# Help AppleClang find Homebrew's libomp if not auto-detected
if(APPLE AND NOT OpenMP_CXX_FOUND)
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



# Make find_package honor UPPER-CASE <PKG>_ROOT variables (CMake 3.30+)
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()
if(POLICY CMP0169)
  cmake_policy(SET CMP0169 NEW)
endif()
if(POLICY CMP0177)
  cmake_policy(SET CMP0177 NEW)
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
endif()

find_package(Threads REQUIRED)
list(APPEND COMMON_LIBRARIES Threads::Threads)


#########################################
# COMPILE BLAS
#########################################
add_subdirectory(dblas)
set (DBLAS_LIBRARIES ${DBLAS_LIBRARIES} dblas)

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

add_library(perclibrary STATIC BaseSpline.cpp MassHandler.cpp PSMDescription.cpp MonotoneRegressor.cpp IsotonicPEP.cpp
								  SetHandler.cpp StdvNormalizer.cpp Caller.cpp CompositionSorter.cpp CrossValidation.cpp Enzyme.cpp Globals.cpp NoNormalizer.cpp Normalizer.cpp
								  SanityCheck.cpp UniNormalizer.cpp DataSet.cpp TabFileValidator.cpp FeatureNames.cpp LogisticRegression.cpp Numerical.cpp Option.cpp PosteriorEstimator.cpp
								  ProteinProbEstimator.cpp Scores.cpp PseudoRandom.cpp Set.cpp ssl.cpp PackedVector.cpp
								  PackedMatrix.cpp MyException.cpp ProteinScoreHolder.cpp PickedProteinInterface.cpp FeatureMemoryPool.cpp Timer.cpp TmpDir.cpp ValidateTabFile.cpp Vector.cpp)

if(OpenMP_CXX_FOUND)
  target_link_libraries(perclibrary PUBLIC ${COMMON_LIBRARIES} OpenMP::OpenMP_CXX)
else()
  target_link_libraries(perclibrary PUBLIC ${COMMON_LIBRARIES})
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
