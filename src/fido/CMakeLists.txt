file(GLOB FIDO_SOURCES Set.cpp Vector.cpp Numerical.cpp Matrix.cpp GramSchmidt.cpp LPSolver.cpp ForwardSampler.cpp Random.cpp ProteinIdentifier.cpp BasicBigraph.cpp BasicGroupBigraph.cpp GroupPowerBigraph.cpp mainFido.cpp PivdoSplitter.cpp)
add_library(fido ${FIDO_SOURCES})

