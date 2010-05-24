link_directories(${PERCOLATOR_BINARY_DIR}/src)

file(GLOB ELUDE_SOURCES *.cpp)

add_executable(elude ${ELUDE_SOURCES})

target_link_libraries(elude perclibrary ArrayLibrary)

install(TARGETS elude EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!
