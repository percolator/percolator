include_directories(${PERCOLATOR_SOURCE_DIR}/src)
link_directories(${PERCOLATOR_SOURCE_DIR}/src)

file(GLOB FIDO_SOURCES Set.cpp Vector.cpp PackedVector.cpp Numerical.cpp Matrix.cpp PackedMatrix.cpp GramSchmidt.cpp LPSolver.cpp ForwardSampler.cpp Random.cpp ProteinIdentifier.cpp BasicBigraph.cpp BasicGroupBigraph.cpp GroupPowerBigraph.cpp mainFido.cpp PivdoSplitter.cpp)
add_library(fido ${FIDO_SOURCES})

# uncomment to add fido to the list of executables to be installed
#add_executable(fido_main mainFido.cpp)
#target_link_libraries(fido_main perclibrary fido)
#install(TARGETS fido_main DESTINATION ./bin/fido)

#add_subdirectory(useful)
