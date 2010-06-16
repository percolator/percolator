link_directories(${PERCOLATOR_BINARY_DIR}/src)

file(GLOB ELUDE_SOURCES *.cpp)

add_executable(elude ${ELUDE_SOURCES})

# we should be able to remove  ${XERCESC_LIBRARIES} ${TOKYOCABINET_LIBRARIES} when we have split out the sqt->xml mzIdentML->xml file format conversion functionality
target_link_libraries(elude perclibrary ArrayLibrary  ${XERCESC_LIBRARIES} ${TOKYOCABINET_LIBRARIES})

install(TARGETS elude EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!
