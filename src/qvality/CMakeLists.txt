link_directories(${PERCOLATOR_BINARY_DIR}/src)

file(GLOB QVALITY_SOURCES *.cpp)

add_executable(qvality ${QVALITY_SOURCES})

target_link_libraries(qvality perclibrary ArrayLibrary)

install(TARGETS qvality EXPORT PERCOLATOR DESTINATION ./bin) # Important to use relative path here (used by CPack)!
