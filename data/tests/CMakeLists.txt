# Defining list with all tests.
set(system_tests_names
  # Correctness tests.
  SystemTest_Percolator_Correctness
  SystemTest_Qvality_Correctness
  # Performance tests.
  SystemTest_Percolator_Performance
  SystemTest_Qvality_Performance
)

IF(NOT EXCLUDE_ELUDE)
  set(system_tests_names "${system_tests_names}, SystemTest_Elude_Performance")
ENDIF()

set(system_tests_dir "${CMAKE_SOURCE_DIR}/data/tests")

# NOTE: To run system tests (after building), run 'make test' in the 'output' folder.
foreach(current_system_test ${system_tests_names})
  #MESSAGE ("Copying test: ${current_system_test} ...")
  add_test(${current_system_test} python ${system_tests_dir}/${current_system_test}.py)
endforeach(current_system_test)


# Google tests are built explicitely.
IF (GOOGLE_TEST)
  #find_package(GTest REQUIRED)
  include_directories (${PERCOLATOR_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/../gtest-1.5.0/include/)
  add_executable (gtest_unit simple_tests.cpp)
  target_link_libraries (gtest_unit perclibrary ArrayLibrary ${XERCESC_LIBRARIES} ${CMAKE_BINARY_DIR}/../gtest-1.5.0/build/libgtest.a -pthread)
  IF (NOT EXCLUDE_CONVERTERS)
    target_link_libraries (gtest_unit perclibrary ArrayLibrary MSToolkit ${XERCESC_LIBRARIES} ${CMAKE_BINARY_DIR}/../gtest-1.5.0/build/libgtest.a -pthread)
  ENDIF ()
  install (TARGETS gtest_unit EXPORT PERCOLATOR DESTINATION .) # Important to use relative path here (used by CPack)!
ENDIF(GOOGLE_TEST)
