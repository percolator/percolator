# Adding CTest basic system-level tests here.
# IF(SYSTEM_TESTS STREQUAL TRUE)

  # Defining list with all tests.
  set(system_tests_names
      # Testing help manuals.
#      SystemTest_Percolator_Help
#      SystemTest_Qvality_Help
#      SystemTest_Elude_Help
      
      # Basic functionality tests.
#      SystemTest_Percolator_BasicFunc
#      SystemTest_Qvality_BasicFunc
#      SystemTest_Elude_BasicFunc
      
      # Correctness tests.
#      SystemTest_Percolator_Correctness
#      SystemTest_Qvality_Correctness
#      SystemTest_Elude_Correctness
      
      # Performance tests.
      SystemTest_Percolator_Performance
#      SystemTest_Qvality_Performance
#      SystemTest_Elude_Performance
      
      # Workload and other types of tests.
#      SystemTest_Percolator_Workload
#      SystemTest_Qvality_Workload
#      SystemTest_Elude_Workload
    )
  # Just copying the script (permissions are preserved).
  set(system_tests_dir "${CMAKE_BINARY_DIR}/../bin")

  # NOTE: To run system tests (after running build_me.sh), run 'make test' in the 'output' folder.
  foreach(current_system_test ${system_tests_names})
    MESSAGE ("Copying test: ${current_system_test} ...")
    configure_file (${current_system_test}.py ${system_tests_dir}/${current_system_test}.py COPYONLY)
    add_test(${current_system_test} python ${system_tests_dir}/${current_system_test}.py)
  endforeach(current_system_test)

# ENDIF(SYSTEM_TESTS STREQUAL TRUE)

# Google tests are built explicitely.
IF (GOOGLE_TEST STREQUAL TRUE)
  include_directories (${PERCOLATOR_SOURCE_DIR}/src)
  add_executable (gtest_unit simple_tests.cpp)
  target_link_libraries (gtest_unit libgtest.a perclibrary ArrayLibrary MSToolkit)
  install (TARGETS gtest_unit EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!
ENDIF (GOOGLE_TEST STREQUAL TRUE)
