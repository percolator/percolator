# SETTING VARIBALES THAT WILL BE OF USE IN THE TEST SCRIPTS
#set(pathToBinaries ${CMAKE_INSTALL_PREFIX}/bin)
#set(pathToData ${CMAKE_SOURCE_DIR}/data)
#set(pathToWritable ${ELUDE_WRITABLE_DIR})
set(pathToBinaries ${CMAKE_BINARY_DIR})
set(pathToData ${CMAKE_SOURCE_DIR}/data)
set(pathToWritable ${CMAKE_BINARY_DIR}/Testing)

# STORE NEWLY SET VARIABLES IN *.h.cmake FILES
file(GLOB_RECURSE configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/${myReplFile} @ONLY)
endforeach()

set(tests
  SystemTest_Elude_TrainEvaluate
  SystemTest_Elude_AutomaticModel
  SystemTest_Elude_SaveLoadEvaluate
  SystemTest_Elude_PtmsTrainEvaluate
  SystemTest_Elude_PtmsAutomaticModel
  SystemTest_Elude_PtmsSaveLoad
)

set(system_tests_names ${tests})
set(system_tests_dir ${CMAKE_CURRENT_BINARY_DIR})

# TO RUN THE TESTS: BUILD, COMPILE, INSTALL AND INVOKE 'make test' FROM BUILD FOLDER
foreach(current_system_test ${system_tests_names})
  add_test(${current_system_test} python ${system_tests_dir}/${current_system_test}.py)
endforeach(current_system_test)
