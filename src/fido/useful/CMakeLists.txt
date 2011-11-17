include_directories(${PERCOLATOR_SOURCE_DIR}/src/fido)
link_directories(${PERCOLATOR_SOURCE_DIR}/src/fido)

add_executable(fdrDivergence fdrDivergence.cpp)
target_link_libraries(fdrDivergence perclibrary fido)
install(TARGETS fdrDivergence DESTINATION ./bin/fido)

add_executable(rlistFDR rlistFDR.cpp)
target_link_libraries(rlistFDR perclibrary fido)
install(TARGETS rlistFDR DESTINATION ./bin/fido)

add_executable(rlistROC rlistROC.cpp)
target_link_libraries(rlistROC perclibrary fido)
install(TARGETS rlistROC DESTINATION ./bin/fido)

add_executable(roc50 roc50.cpp)
target_link_libraries(roc50 perclibrary fido)
install(TARGETS roc50 DESTINATION ./bin/fido)
