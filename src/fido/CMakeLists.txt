include_directories(${PERCOLATOR_SOURCE_DIR}/src)
link_directories(${PERCOLATOR_SOURCE_DIR}/src)

file(GLOB FIDO_SOURCES Set.cpp Vector.cpp Numerical.cpp Random.cpp BasicBigraph.cpp BasicGroupBigraph.cpp GroupPowerBigraph.cpp)
add_library(fido ${FIDO_SOURCES})

