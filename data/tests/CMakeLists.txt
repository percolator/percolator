set(tests
  SystemTest_Percolator_Correctness
  SystemTest_Qvality_Correctness
  SystemTest_Percolator_Performance
  SystemTest_Qvality_Performance
  UnitTest_Percolator_RunAllTests
)
# Defining list with all tests.
#IF(EXCLUDE_ELUDE AND EXCLUDE_CONVERTERS)
#set(tests
#  SystemTest_Percolator_Correctness
#  SystemTest_Qvality_Correctness
#  SystemTest_Percolator_Performance
#  SystemTest_Qvality_Performance
#  UnitTest_Percolator_RunAllTests
#)
#ELSEIF((NOT EXCLUDE_ELUDE) AND EXCLUDE_CONVERTERS)
#set(tests
#  SystemTest_Percolator_Correctness
#  SystemTest_Qvality_Correctness
#  SystemTest_Percolator_Performance
#  SystemTest_Qvality_Performance
#  SystemTest_Elude_AutomaticModel
#  SystemTest_Elude_PtmsAutomaticModel
#  SystemTest_Elude_PtmsSaveLoad
#  SystemTest_Elude_PtmsTrainEvaluate
#  SystemTest_Elude_SaveLoadEvaluate
#  SystemTest_Elude_TrainEvaluate
#  UnitTest_Elude_RunAllTests
#  UnitTest_Percolator_RunAllTests
#)
#ELSEIF(EXCLUDE_ELUDE AND (NOT EXCLUDE_CONVERTERS))
#set(tests
#  SystemTest_Percolator_Correctness
#  SystemTest_Qvality_Correctness
#  SystemTest_Percolator_Performance
#  SystemTest_Qvality_Performance
#  SystemTest_Sqt2pin_Correctness
#  UnitTest_Percolator_RunAllTests
#)
#ELSE()
#set(tests
#  SystemTest_Sqt2pin_Correctness
#  SystemTest_Percolator_Correctness
#  SystemTest_Qvality_Correctness
#  SystemTest_Percolator_Performance
#  SystemTest_Qvality_Performance
#  SystemTest_Elude_AutomaticModel
#  SystemTest_Elude_PtmsAutomaticModel
#  SystemTest_Elude_PtmsSaveLoad
#  SystemTest_Elude_PtmsTrainEvaluate
#  SystemTest_Elude_SaveLoadEvaluate
#  SystemTest_Elude_TrainEvaluate
#  UnitTest_Elude_RunAllTests
#  UnitTest_Percolator_RunAllTests
#)
#ENDIF()


set(system_tests_names ${tests})
set(system_tests_dir "${CMAKE_SOURCE_DIR}/data/tests")

# NOTE: To run system tests (after building), run 'make test' in the 'output' folder.
foreach(current_system_test ${system_tests_names})
  # MESSAGE ("Copying test: ${current_system_test} ...")
  add_test(${current_system_test} python ${system_tests_dir}/${current_system_test}.py)
endforeach(current_system_test)


# Google tests are built explicitely
  # find_package(GTest REQUIRED)
  include_directories (${PERCOLATOR_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/../../src ${GOOGLE_TEST_PATH}/include/ ${PERCOLATOR_SOURCE_DIR}/src/fido/ ${PERCOLATOR_SOURCE_DIR}/data/tests/)
  add_executable (gtest_unit Unit_tests_Percolator_main.cpp)
  target_link_libraries (gtest_unit perclibrary ${GOOGLE_TEST_PATH}/build/libgtest.a -pthread)
  install (TARGETS gtest_unit EXPORT PERCOLATOR DESTINATION ./bin) # Important to use relative path here (used by CPack)!
#IF (GOOGLE_TEST AND (NOT EXCLUDE_ELUDE))
#  include_directories (${PERCOLATOR_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/../../src ${GOOGLE_TEST_PATH}/include/ ${PERCOLATOR_SOURCE_DIR}/src/fido/ ${PERCOLATOR_SOURCE_DIR}/src/elude/ ${PERCOLATOR_SOURCE_DIR}/data/tests/)
#  add_executable (gtest_unit Unit_tests_Elude_main.cpp)
#  target_link_libraries (gtest_unit eludelibrary ${GOOGLE_TEST_PATH}/build/libgtest.a -pthread)
#  install (TARGETS gtest_unit EXPORT PERCOLATOR DESTINATION ./bin) # Important to use relative path here (used by CPack)!
#ENDIF()
