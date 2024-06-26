# If the user builds without a dependency on GTest, then build nothing.
if(NOT WITHOUT_GTEST)
  # If GTest is not installed, build and install a private copy.
  find_package(GTest)
  if(GTEST_FOUND)
    message(STATUS "Using installed GoogleTest at: ${GTEST_INCLUDE_DIRS}")
  else(GTEST_FOUND)
    message(STATUS "GoogleTest package not found; proceeding to download it.")
    configure_file(GTestExtProject.cmake.in googletest-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
    if(result)
      message(FATAL_ERROR "Failed to download GoogleTest: ${result}.")
    endif()
    message(STATUS "Building GoogleTest in ${CMAKE_CURRENT_BINARY_DIR}/googletest-download.")
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
    if(result)
      message(FATAL_ERROR "Failed to build GoogleTest package: ${result}.")
    endif()
    # Prevent overriding the parent project's compiler/linker
    # settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                     ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                     EXCLUDE_FROM_ALL)
    set(GTEST_INCLUDE_DIRS "${gtest_SOURCE_DIR}/include")
#    set(GTEST_BOTH_LIBRARIES gtest gtestmain)
    set(GTEST_BOTH_LIBRARIES GTest::gtest GTest::gtest_main)
    message(STATUS "GoogleTest has been installed succesfully: ${GTEST_INCLUDE_DIRS}.")
  endif(GTEST_FOUND)

  # Linking and building unit tests
  find_package(Boost ${BOOST_MIN_VERSION} REQUIRED)
  add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
  add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)
  if(WIN32)
    add_definitions(-DBOOST_ALL_NO_LIB) # disable autolinking in boost
  endif(WIN32)
  include_directories(${GTEST_INCLUDE_DIRS}
      ${Boost_INCLUDE_DIRS}
      ${PERCOLATOR_SOURCE_DIR}/src
      # ${PERCOLATOR_SOURCE_DIR}/src/fido
      ${CMAKE_BINARY_DIR}/src)
  add_executable(gtest_unit
      Unit_tests_Percolator_main.cpp
      UnitTest_Percolator_Option.cpp
      UnitTest_Percolator_TabReader.cpp
      UnitTest_Percolator_SetHandler.cpp
      UnitTest_Percolator_DataSet.cpp
      UnitTest_Percolator_Scores.cpp
      UnitTest_Percolator_CrossValidation.cpp)
  # Link with all required libraries
  set(UNIT_TEST_LIBRARIES perclibrary blas)
  if(NOT MSVC)
    set(UNIT_TEST_LIBRARIES ${UNIT_TEST_LIBRARIES} pthread)
    if(APPLE)
      set(UNIT_TEST_LIBRARIES ${UNIT_TEST_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
    endif(APPLE)
  endif(NOT MSVC)
  if(GTEST_FOUND)
    set(UNIT_TEST_LIBRARIES ${UNIT_TEST_LIBRARIES} ${GTEST_BOTH_LIBRARIES})
  else(GTEST_FOUND)
    set(UNIT_TEST_LIBRARIES ${UNIT_TEST_LIBRARIES} gtest gtest_main)
  endif(GTEST_FOUND)
  if(COVERAGE)
    # Insert flags for generating coverage data
    target_compile_options(perclibrary PUBLIC -ftest-coverage -fprofile-arcs)
    target_compile_options(gtest_unit PUBLIC -ftest-coverage -fprofile-arcs)
    set(UNIT_TEST_LIBRARIES -fprofile-arcs ${UNIT_TEST_LIBRARIES})
  endif(COVERAGE)
  target_link_libraries(gtest_unit ${UNIT_TEST_LIBRARIES})
  add_test(UnitTest_Percolator_RunAllTests gtest_unit)
  # Important to use relative paths here (used by CPack)!
  install(TARGETS gtest_unit EXPORT PERCOLATOR DESTINATION ./bin)
endif(NOT WITHOUT_GTEST)
