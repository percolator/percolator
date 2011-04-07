# SETTING VARIBALES THAT WILL BE OF USE IN THE TEST SCRIPTS
set(pathToBinaries ${CMAKE_INSTALL_PREFIX}/bin)

# STORE NEWLY SET VARIABLES IN *.h.cmake FILES
file(GLOB_RECURSE configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/${myReplFile} @ONLY)
endforeach()

set(system_tests_names UnitTest_Percolator_RunAllTests)
set(system_tests_dir ${CMAKE_CURRENT_BINARY_DIR})

# TO RUN THE TESTS: BUILD, COMPILE, INSTALL AND INVOKE 'make test' FROM BUILD FOLDER
foreach(current_system_test ${system_tests_names})
  add_test(${current_system_test} python ${system_tests_dir}/${current_system_test}.py)
endforeach(current_system_test)

# USE THIS WHEN THE 1.5 VERSION OF GTEST WILL BE FINALLY SUPPORTED INSTEAD OF BUILDING IN INSTALL SCRIPT
# find_package(GTest REQUIRED)

# LINING AND BUILDING GTEST
include_directories (${GOOGLE_TEST_PATH}/include ${PERCOLATOR_SOURCE_DIR}/src ${PERCOLATOR_SOURCE_DIR}/src/fido ${PERCOLATOR_SOURCE_DIR}/data/tests ${CMAKE_BINARY_DIR}/src)
add_executable (gtest_unit Unit_tests_Percolator_main.cpp)
target_link_libraries (gtest_unit perclibrary ${GOOGLE_TEST_PATH}/build/libgtest.a -pthread)
install (TARGETS gtest_unit EXPORT PERCOLATOR DESTINATION ./bin) # Important to use relative path here (used by CPack)!

#IF (GOOGLE_TEST AND (NOT EXCLUDE_ELUDE))
#  include_directories (${PERCOLATOR_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/../../src ${GOOGLE_TEST_PATH}/include/ ${PERCOLATOR_SOURCE_DIR}/src/fido/ ${PERCOLATOR_SOURCE_DIR}/src/elude/ ${PERCOLATOR_SOURCE_DIR}/data/tests/)
#  add_executable (gtest_unit Unit_tests_Elude_main.cpp)
#  target_link_libraries (gtest_unit eludelibrary ${GOOGLE_TEST_PATH}/build/libgtest.a -pthread)
#  install (TARGETS gtest_unit EXPORT PERCOLATOR DESTINATION ./bin) # Important to use relative path here (used by CPack)!
#ENDIF()
