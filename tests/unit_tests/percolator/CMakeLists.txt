set(CMAKE_CXX_STANDARD 14)

include(FetchContent)

if(WIN32)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()

set(BUILD_SHARED_LIBS OFF)

# =============================
# Fetch or find GoogleTest
# =============================
message(STATUS "Fetching GoogleTest...")
set(BUILD_GTEST ON)
set(INSTALL_GTEST OFF)
set(BUILD_GMOCK OFF)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.16.0
)
FetchContent_Populate(googletest)

# This is the correct way to build the actual library
add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})

# Optional: apply RTTI to gtest explicitly (not usually needed if add_compile_options is set globally)
# target_compile_options(gtest PRIVATE -frtti)
# target_compile_options(gtest_main PRIVATE -frtti)

enable_testing()
include(GoogleTest)

# =============================
# Fetch or find Eigen3
# =============================
find_package(Eigen3 QUIET)
if(NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found, downloading it...")
    FetchContent_Declare(
      eigen
      GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
      GIT_TAG 3.4.0
    )
    FetchContent_MakeAvailable(eigen)
endif()

# =============================
# Find Boost
# =============================
find_package(Boost ${BOOST_MIN_VERSION} REQUIRED)
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)
if(WIN32)
    add_definitions(-DBOOST_ALL_NO_LIB)
endif()

# =============================
# Add test executable
# =============================
add_executable(gtest_unit
    Unit_tests_Percolator_main.cpp
    UnitTest_Percolator_Option.cpp
    UnitTest_Percolator_TabReader.cpp
    UnitTest_Percolator_SetHandler.cpp
    UnitTest_Percolator_DataSet.cpp
    UnitTest_Percolator_IsplineRegression.cpp
    UnitTest_Percolator_Scores.cpp
    UnitTest_Percolator_CrossValidation.cpp
)

# =============================
# Include directories
# =============================
target_include_directories(gtest_unit
  PRIVATE
    ${Boost_INCLUDE_DIRS}
    ${PERCOLATOR_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/src
)

# =============================
# Link libraries
# =============================
target_link_libraries(gtest_unit
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    Eigen3::Eigen
    perclibrary
    dblas
)
if(NOT WIN32)
  target_link_libraries(gtest_unit PRIVATE pthread)
endif()

# (Optional: if you want coverage)
if(COVERAGE)
    target_compile_options(gtest_unit PRIVATE -ftest-coverage -fprofile-arcs)
    target_link_options(gtest_unit PRIVATE -ftest-coverage -fprofile-arcs)
endif()

# =============================
# Enable test
# =============================
add_test(NAME UnitTest_Percolator_RunAllTests COMMAND gtest_unit)

# =============================
# Install
# =============================
install(TARGETS gtest_unit EXPORT PERCOLATOR DESTINATION ./bin)
