# SETTING VARIBALES THAT WILL BE OF USE IN THE TEST SCRIPTS
set(pathToBinaries ${CMAKE_INSTALL_PREFIX}/bin)
set(pathToData ${CMAKE_SOURCE_DIR}/data)
set(pathToOutputData ${CMAKE_BINARY_DIR}/data)
if(XML_SUPPORT)
  set(xmlSupport True)
else(XML_SUPPORT)
  set(xmlSupport False)
endif(XML_SUPPORT)

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
              COMMAND python ${system_tests_dir}/${current_system_test}.py)
endforeach(current_system_test)
add_custom_target(test-install ${system_tests_commands} VERBATIM)
