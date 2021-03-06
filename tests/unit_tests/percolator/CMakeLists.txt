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
    set(GTEST_BOTH_LIBRARIES gtest gtestmain)
    message(STATUS "GoogleTest has been installed succesfully: ${GTEST_INCLUDE_DIRS}.")
  endif(GTEST_FOUND)

  # Linking and building unit tests
  include_directories(${GTEST_INCLUDE_DIRS}
      ${PERCOLATOR_SOURCE_DIR}/src
      ${PERCOLATOR_SOURCE_DIR}/src/fido
      ${CMAKE_BINARY_DIR}/src)
  add_executable(gtest_unit
      Unit_tests_Percolator_main.cpp
      UnitTest_Percolator_Fido.cpp
      UnitTest_Percolator_Option.cpp
      UnitTest_Percolator_TabReader.cpp
      UnitTest_Percolator_SetHandler.cpp
      UnitTest_Percolator_DataSet.cpp
      UnitTest_Percolator_Scores.cpp)
  # Flags for generating coverage data
  if(COVERAGE)
    target_compile_options(perclibrary PUBLIC -ftest-coverage -fprofile-arcs)
    target_compile_options(gtest_unit PUBLIC -ftest-coverage -fprofile-arcs)
    target_link_libraries(gtest_unit -fprofile-arcs)
  endif(COVERAGE)
  if (GTEST_FOUND)
    target_link_libraries(gtest_unit perclibrary fido pthread ${GTEST_BOTH_LIBRARIES})
  else (GTEST_FOUND)
    target_link_libraries(gtest_unit perclibrary fido pthread gtest gtest_main)
  endif(GTEST_FOUND)
  add_test(UnitTest_Percolator_RunAllTests gtest_unit)

  # Important to use relative paths here (used by CPack)!
  install(TARGETS gtest_unit EXPORT PERCOLATOR DESTINATION ./bin)
endif(WITHOUT_GTEST)
