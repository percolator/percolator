link_directories(${PERCOLATOR_BINARY_DIR}/src)

file(GLOB QVALITY_SOURCES *.cpp)

add_executable(qvality ${QVALITY_SOURCES})

if(COVERAGE)
  target_link_libraries(qvality -fprofile-arcs)
endif(COVERAGE)
target_link_libraries(qvality perclibrary)

install(TARGETS qvality EXPORT PERCOLATOR DESTINATION ./bin) # Important to use relative path here (used by CPack)!
