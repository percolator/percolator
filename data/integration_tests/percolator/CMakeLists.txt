
set(pathToTestScripts ${CMAKE_BINARY_DIR}/data/integration_tests/percolator)
set(pathToOutputData ${CMAKE_BINARY_DIR}/data/integration_tests/percolator/testData)
set(pathToBinaries ${CMAKE_BINARY_DIR}/src)


file(MAKE_DIRECTORY ${pathToOutputData})

file(GLOB_RECURSE configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.py )
foreach(i ${configurefiles})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${i} ${CMAKE_CURRENT_BINARY_DIR}/${i} @ONLY)
endforeach()

