add_executable(qvality qvality_main.cpp)

if(COVERAGE)
  target_link_libraries(qvality -fprofile-arcs)
endif(COVERAGE)
target_link_libraries(qvality perclibrary)

install(TARGETS qvality EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!
