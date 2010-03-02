link_directories(${LUMINATOR_BINARY_DIR}/src)

set(QVALITY_ALIAS "qvalumity")

file(GLOB QVALITY_SOURCES *.cpp)

add_executable(${QVALITY_ALIAS} ${QVALITY_SOURCES})

target_link_libraries(${QVALITY_ALIAS} perclibrary ArrayLibrary)

install(TARGETS ${QVALITY_ALIAS} EXPORT ${PROJECT_NAME} DESTINATION bin) # Important to use relative path here (used by CPack)!
