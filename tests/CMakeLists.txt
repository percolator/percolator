include_directories (${PERCOLATOR_SOURCE_DIR}/src)

add_executable (gtest_unit simple_tests.cpp)

target_link_libraries (gtest_unit libgtest.a perclibrary ArrayLibrary MSToolkit)

install (TARGETS gtest_unit EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!

add_test (SimplePercTest1 SimpleTest1.sh "--help")
add_test (SimpleQvalTest2 SimpleTest2.sh "--help")

configure_file(SimpleTest1.sh ${CMAKE_BINARY_DIR}/tests/SimpleTest1.sh COPYONLY)
configure_file(SimpleTest2.sh ${CMAKE_BINARY_DIR}/tests/SimpleTest2.sh COPYONLY)