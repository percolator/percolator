set(CMAKE_CXX_STANDARD 17)

# =============================
# Find or fetch GoogleTest
# =============================
# An installed GoogleTest (e.g. libgtest-dev or gtest-devel) is used when one
# is found, so that configuring needs no network access. Otherwise GoogleTest
# is downloaded with FetchContent. On MSVC an installed GoogleTest may use
# another C runtime than percolator, so it is used only on request there.
if(MSVC)
  set(_use_system_gtest_default OFF)
else()
  set(_use_system_gtest_default ON)
endif()
option(USE_SYSTEM_GTEST "Use an installed GoogleTest if one is found"
       ${_use_system_gtest_default})
if(USE_SYSTEM_GTEST)
  find_package(GTest CONFIG QUIET)
endif()

if(GTest_FOUND)
  message(STATUS "Using GoogleTest ${GTest_VERSION} from ${GTest_DIR}")
else()
  message(STATUS "Fetching GoogleTest...")
  include(FetchContent)

  if(WIN32)
    # GoogleTest then follows CMAKE_MSVC_RUNTIME_LIBRARY (set in the top
    # CMakeLists.txt) instead of rewriting the compiler flags itself
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  endif()

  set(BUILD_SHARED_LIBS OFF)
  set(BUILD_GTEST ON)
  set(INSTALL_GTEST OFF)
  set(BUILD_GMOCK OFF)

  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.16.0
  )
  FetchContent_MakeAvailable(googletest)
endif()

# Optional: apply RTTI to gtest explicitly (not usually needed if add_compile_options is set globally)
# target_compile_options(gtest PRIVATE -frtti)
# target_compile_options(gtest_main PRIVATE -frtti)

enable_testing()
include(GoogleTest)

# =============================
# Add test executable
# =============================
add_executable(gtest_unit
    Unit_tests_Percolator_main.cpp
    UnitTest_Percolator_BlasLevel1.cpp
    UnitTest_Percolator_Option.cpp
    UnitTest_Percolator_TabReader.cpp
    UnitTest_Percolator_SetHandler.cpp
    UnitTest_Percolator_DataSet.cpp
    UnitTest_Percolator_IsplineRegression.cpp
    UnitTest_Percolator_InferPep.cpp
    UnitTest_Percolator_Scores.cpp
    UnitTest_Percolator_CrossValidation.cpp
    UnitTest_Percolator_Normalizer.cpp
    UnitTest_Percolator_PSMDescription.cpp
    UnitTest_Percolator_PosteriorEstimator.cpp
    UnitTest_Percolator_SanityCheck.cpp
    UnitTest_Percolator_Caller.cpp
    UnitTest_Percolator_CompositionSorter.cpp
    UnitTest_Percolator_TrainingMemory.cpp
    UnitTest_Percolator_ParallelSort.cpp
    UnitTest_Percolator_ProteinDigest.cpp
)

# =============================
# Include directories
# =============================
target_include_directories(gtest_unit
  PRIVATE
    ${PERCOLATOR_SOURCE_DIR}/src
    ${PERCOLATOR_SOURCE_DIR}/src/picked_protein
    ${CMAKE_BINARY_DIR}/src
)

# =============================
# Link libraries
# =============================
target_link_libraries(gtest_unit
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    perclibrary
    picked_protein
)
# GoogleTest needs the threads library; percolator itself uses OpenMP only
find_package(Threads REQUIRED)
target_link_libraries(gtest_unit PRIVATE Threads::Threads)

# (Optional: if you want coverage)
if(COVERAGE)
    target_compile_options(gtest_unit PRIVATE -ftest-coverage -fprofile-arcs)
    target_link_options(gtest_unit PRIVATE -ftest-coverage -fprofile-arcs)
endif()

# =============================
# Enable test
# =============================
add_test(NAME UnitTest_Percolator_RunAllTests COMMAND gtest_unit)
