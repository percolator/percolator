

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
  # Homebrew prefixes on Apple silicon and on Intel
  foreach(_HOMEBREW_OMP_ROOT "/opt/homebrew/opt/libomp" "/usr/local/opt/libomp")
    if (EXISTS "${_HOMEBREW_OMP_ROOT}/include" AND EXISTS "${_HOMEBREW_OMP_ROOT}/lib")
      # Link the static library where there is one, so that the binary
      # runs without Homebrew
      set(_HOMEBREW_OMP_LIBRARY "${_HOMEBREW_OMP_ROOT}/lib/libomp.a")
      if (NOT EXISTS "${_HOMEBREW_OMP_LIBRARY}")
        set(_HOMEBREW_OMP_LIBRARY "${_HOMEBREW_OMP_ROOT}/lib/libomp.dylib")
      endif()
      # Hint include/lib to CMake's OpenMP finder
      set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp" CACHE STRING "" FORCE)
      set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "" FORCE)
      set(OpenMP_omp_LIBRARY "${_HOMEBREW_OMP_LIBRARY}" CACHE FILEPATH "" FORCE)
      include_directories("${_HOMEBREW_OMP_ROOT}/include")
      # Retry discovery
      find_package(OpenMP)
      break()
    endif()
  endforeach()
endif()
if (OpenMP_CXX_FOUND)
    list(APPEND COMMON_LIBRARIES OpenMP::OpenMP_CXX)
    # With STATIC on Linux, FindOpenMP takes libgomp.a, which calls dlopen
    # and dlsym. Before glibc 2.34 they are in libdl, not in libc (the
    # manylinux_2_28 wheels and the portable build), so libdl follows it.
    if(UNIX AND NOT APPLE AND STATIC)
      list(APPEND COMMON_LIBRARIES ${CMAKE_DL_LIBS})
    endif()
else()
  message(WARNING "OpenMP was not found, so percolator will run on one "
                  "thread only. On macOS, install libomp (brew install libomp).")
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


if (MSVC)
  add_compile_definitions(
    _HAS_STD_BYTE=0     # disables std::byte so 'byte' in Windows headers is unambiguous
    NOMINMAX            # avoid macro min/max conflicts
    _WIN32_WINNT=0x0601 # silence that SDK warning (Windows 7 or higher)
  )
endif()

#########################################
# STANDARD LIBRARY FEATURES
#########################################
# std::to_chars and std::from_chars for double are missing from older
# standard libraries (GCC before 11), and libc++ provides them on macOS only
# from some deployment targets on, without announcing them with
# __cpp_lib_to_chars. The code falls back to slower equivalents without them.
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
#include <charconv>
int main() {
  char buf[32];
  double value = 0.5;
  return std::to_chars(buf, buf + sizeof(buf), value,
                       std::chars_format::general, 6).ec == std::errc() ? 0 : 1;
}" PERCOLATOR_HAVE_TO_CHARS_DOUBLE)
check_cxx_source_compiles("
#include <charconv>
int main() {
  const char text[] = \"0.5\";
  double value = 0.0;
  return std::from_chars(text, text + 3, value).ec == std::errc() ? 0 : 1;
}" PERCOLATOR_HAVE_FROM_CHARS_DOUBLE)
set(PERCOLATOR_CHARCONV_DEFINITIONS "")
foreach(conversion TO_CHARS FROM_CHARS)
  if(PERCOLATOR_HAVE_${conversion}_DOUBLE)
    list(APPEND PERCOLATOR_CHARCONV_DEFINITIONS PERCOLATOR_${conversion}_DOUBLE=1)
  else()
    list(APPEND PERCOLATOR_CHARCONV_DEFINITIONS PERCOLATOR_${conversion}_DOUBLE=0)
  endif()
endforeach()
# for every target in this directory and below; perclibrary passes them on
# to the targets that link it, such as the unit tests
add_compile_definitions(${PERCOLATOR_CHARCONV_DEFINITIONS})

#############################################################################
# COMPILE PERCOLATOR
###############################################################################

add_library(perclibrary STATIC BaseSpline.cpp MassHandler.cpp PSMDescription.cpp MonotoneRegressor.cpp IsotonicPEP.cpp
								  SetHandler.cpp StdvNormalizer.cpp Caller.cpp CompositionSorter.cpp CrossValidation.cpp Enzyme.cpp Globals.cpp NoNormalizer.cpp Normalizer.cpp
								  SanityCheck.cpp UniNormalizer.cpp DataSet.cpp TabFileValidator.cpp FeatureNames.cpp LogisticRegression.cpp Numerical.cpp Option.cpp PosteriorEstimator.cpp
								  ProteinProbEstimator.cpp Scores.cpp PseudoRandom.cpp Set.cpp ssl.cpp PackedVector.cpp
								  PackedMatrix.cpp MyException.cpp ProteinScoreHolder.cpp PickedProteinInterface.cpp FeatureMemoryPool.cpp Timer.cpp TmpDir.cpp ValidateTabFile.cpp Vector.cpp
								  picked_protein/PickedProteinCaller.cpp picked_protein/ProteinDigest.cpp)

if(OpenMP_CXX_FOUND)
  target_link_libraries(perclibrary PUBLIC ${COMMON_LIBRARIES} OpenMP::OpenMP_CXX)
else()
  target_link_libraries(perclibrary PUBLIC ${COMMON_LIBRARIES})
endif()
target_compile_definitions(perclibrary INTERFACE ${PERCOLATOR_CHARCONV_DEFINITIONS})

# The picked-protein sources (src/picked_protein) are part of perclibrary:
# they and PickedProteinInterface.cpp use each other, so two static
# libraries had to be linked in one order. The target picked_protein remains
# as a name for targets that link "perclibrary picked_protein"; it builds no
# library file.
add_library(picked_protein INTERFACE)
target_link_libraries(picked_protein INTERFACE perclibrary)

#INCLUDE PERCOLATOR HEADERS
include_directories(.)

#INCLUDE PICKED PROTEIN HEADERS
include_directories(picked_protein)

################################################################################

# COMPILE PERCOLATOR
if(CRUX)
  ADD_DEFINITIONS(-DCRUX)
  include_directories(${CRUX} ${CRUX}/src ${EXT_BINARY_DIR})
ELSE(CRUX)
  add_executable(percolator main.cpp)

  # SET LIBRARIES FOR PERCOLATOR
  list(APPEND COMMON_LIBRARIES perclibrary picked_protein)

  target_link_libraries(percolator PRIVATE ${COMMON_LIBRARIES})

  # INSTALL PERCOLATOR
  install(TARGETS percolator EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!
ENDIF(CRUX)

###############################################################################
# COMPILE QVALITY
###############################################################################

add_subdirectory(qvality)
