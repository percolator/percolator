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
SET( BASE "/tmp")
if(GOOGLE_TEST) ## DO IT IF ONLY IF GOOGLE TEST IS ACTIVATED
  find_package(GoogleTest REQUIRED) 
  if(NOT GTEST_FOUND)
    SET( GT_SITE "http://googletest.googlecode.com/files")
    SET( GT_PACK "gtest-1.6.0")
    message( STATUS "Package Google Test has not been found, prodecing to download it")

    execute_process(COMMAND wget -nc "${GT_SITE}/${GT_PACK}.zip" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
    execute_process(COMMAND unzip -o "${GT_PACK}.zip" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
    execute_process(COMMAND rm "${GT_PACK}.zip" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)

    if(${rv} EQUAL 0) #installed okay
      message(STATUS "GoogleTest has been downloaded succesfully : ${BASE}/${GT_PACK}")
    else(${rv} EQUAL 0)
      message(FATAL_ERROR "There has been a problem downloading the package.Error: '${error}")
    endif(${rv} EQUAL 0)

    message(STATUS "Installing GoogleTest")

    execute_process(COMMAND cmake . WORKING_DIRECTORY "${BASE}/${GT_PACK}" OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
    execute_process(COMMAND make WORKING_DIRECTORY "${BASE}/${GT_PACK}" OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
    #TODO compile it properly and delete it after compiling
    if(${rv} EQUAL 0) #compiled okay
      set(GOOGLE_TEST_PATH "${BASE}/${GT_PACK}") ##set the variable accordingly
      find_package(GoogleTest REQUIRED)
      message(STATUS "GoogleTest has been installed succesfully ${GTEST_INCLUDE_DIRS}")
      
    else(${rv} EQUAL 0)
      message(FATAL_ERROR "The package Google Test has not been found, you can download it and install it from :\n
			  http://code.google.com/p/googletest/downloads/list .Alternatively you can run this instruction to install it\n
			  under debian based systems : apt-get install google-gadgets-gtk (you must run it as sudo).")
    endif(${rv} EQUAL 0)

  else(NOT GTEST_FOUND)
    message(STATUS "Package Google Test found: ${GTEST_INCLUDE_DIRS}")
  endif(NOT GTEST_FOUND)

endif(GOOGLE_TEST)

# LINING AND BUILDING GTEST
include_directories (${GTEST_INCLUDE_DIRS} ${PERCOLATOR_SOURCE_DIR}/src ${PERCOLATOR_SOURCE_DIR}/src/fido ${PERCOLATOR_SOURCE_DIR}/data/tests ${CMAKE_BINARY_DIR}/src)
add_executable (gtest_unit Unit_tests_Percolator_main.cpp)
target_link_libraries (gtest_unit perclibrary ${GTEST_BOTH_LIBRARIES} pthread)
add_test(AllTestsInFoo gtest_unit)
install (TARGETS gtest_unit EXPORT PERCOLATOR DESTINATION ./bin) # Important to use relative path here (used by CPack)!
