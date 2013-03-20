#include_directories(${PERCOLATOR_SOURCE_DIR}/src)
#link_directories(${PERCOLATOR_SOURCE_DIR}/src)
#link_directories(${PERCOLATOR_BINARY_DIR}/src)

#file(GLOB FIDO_SOURCES Set.cpp Vector.cpp Numerical.cpp Random.cpp BasicBigraph.cpp BasicGroupBigraph.cpp GroupPowerBigraph.cpp)

#add_executable(Fido ${FIDO_SOURCES})

#target_link_libraries(Fido perclibrary)

#install(TARGETS Fido EXPORT PERCOLATOR DESTINATION ./bin) # Important to use relative path here (used by CPack)!

#add_library(fido STATIC ${FIDO_SOURCES})

include_directories(${PERCOLATOR_SOURCE_DIR}/src)
link_directories(${PERCOLATOR_SOURCE_DIR}/src)

file(GLOB FIDO_SOURCES Set.cpp Vector.cpp Numerical.cpp Random.cpp BasicBigraph.cpp BasicGroupBigraph.cpp GroupPowerBigraph.cpp)
#add_library(fido ${FIDO_SOURCES})
add_library(fido STATIC ${FIDO_SOURCES})
