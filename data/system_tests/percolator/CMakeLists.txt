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
# Defining list with all tests.
#IF(EXCLUDE_ELUDE AND EXCLUDE_CONVERTERS)
#set(tests
#  SystemTest_Percolator_Correctness
#  SystemTest_Qvality_Correctness
#  SystemTest_Percolator_Performance
#  SystemTest_Qvality_Performance
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
#)
#ELSEIF(EXCLUDE_ELUDE AND (NOT EXCLUDE_CONVERTERS))
#set(tests
#  SystemTest_Percolator_Correctness
#  SystemTest_Qvality_Correctness
#  SystemTest_Percolator_Performance
#  SystemTest_Qvality_Performance
#  SystemTest_Sqt2pin_Correctness
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
#)
#ENDIF()

set(system_tests_names ${tests})
set(system_tests_dir ${CMAKE_CURRENT_BINARY_DIR})

# Temporarily disabled so that "make test" only runs unit tests.
# TODO(breadbox): Give this section its own build target.
#
#foreach(current_system_test ${system_tests_names})
#  add_test(${current_system_test} python ${system_tests_dir}/${current_system_test}.py)
#endforeach(current_system_test)
