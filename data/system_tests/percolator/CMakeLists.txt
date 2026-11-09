# Find a suitable Python interpreter; without one there are no system tests
find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
  message(STATUS "Python 3 not found: the system tests are not available")
  return()
endif()

# SETTING VARIABLES THAT WILL BE OF USE IN THE TEST SCRIPTS
set(pathToBinaries ${CMAKE_INSTALL_PREFIX}/bin)
set(pathToData ${CMAKE_SOURCE_DIR}/data)
set(pathToOutputData ${CMAKE_BINARY_DIR}/data)
set(xmlSupport False)

# STORE NEWLY SET VARIABLES IN *.h.cmake FILES
file(GLOB_RECURSE configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/${myReplFile} @ONLY)
endforeach()

set(tests
  SystemTest_Percolator_Correctness
  SystemTest_Percolator_Performance
  SystemTest_Qvality_Correctness
  SystemTest_Qvality_Performance
)

# Define the "test-install" target to run the system test scripts
# against the percolator installation.
set(system_tests_dir ${CMAKE_CURRENT_BINARY_DIR})
set(system_tests_commands)
foreach(current_system_test ${tests})
  list(APPEND system_tests_commands
              COMMAND  ${Python3_EXECUTABLE} ${system_tests_dir}/${current_system_test}.py)
endforeach(current_system_test)
add_custom_target(test-install ${system_tests_commands} VERBATIM)

# The same scripts also run with ctest, on the binaries in the build tree,
# so that no installation (and no sudo) is needed: SYSTEM_TEST_BIN_DIR
# replaces pathToBinaries. The performance tests read the output of the
# correctness tests, so they run after them.
foreach(current_system_test ${tests})
  if(current_system_test MATCHES "Qvality")
    set(binary qvality)
  else()
    set(binary percolator)
  endif()
  if(TARGET ${binary})
    add_test(NAME ${current_system_test}
             COMMAND ${CMAKE_COMMAND} -E env
                     "SYSTEM_TEST_BIN_DIR=$<TARGET_FILE_DIR:${binary}>"
                     ${Python3_EXECUTABLE}
                     ${system_tests_dir}/${current_system_test}.py)
  endif()
endforeach()
if(TARGET percolator)
  set_tests_properties(SystemTest_Percolator_Performance PROPERTIES
                       DEPENDS SystemTest_Percolator_Correctness)
endif()
set_tests_properties(SystemTest_Qvality_Performance PROPERTIES
                     DEPENDS SystemTest_Qvality_Correctness)
