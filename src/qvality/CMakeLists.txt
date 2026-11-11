add_executable(qvality qvality_main.cpp)

target_link_libraries(qvality perclibrary)

install(TARGETS qvality EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!
