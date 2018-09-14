include_directories(${PERCOLATOR_SOURCE_DIR}/src)
link_directories(${PERCOLATOR_SOURCE_DIR}/src)

set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> rcv <TARGET> <LINK_FLAGS> <OBJECTS>")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -Wconversion -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wconversion -fPIC")

file(GLOB BLAS_SOURCES dscal.c daxpy.c ddot.c dnrm2.c dgemv.c)
add_library(blas STATIC ${BLAS_SOURCES})
