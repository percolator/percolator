


file(GLOB_RECURSE configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.py )
foreach(i ${configurefiles})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${i} ${CMAKE_CURRENT_BINARY_DIR}/${i} @ONLY)
endforeach()