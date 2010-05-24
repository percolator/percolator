set(CMAKE_CXX_FLAGS "-O3 -static")
set(CMAKE_C_FLAGS "-O3 -static")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -static")
set(CMAKE_C_FLAGS_RELEASE "-O3 -static")

set(CMAKE_CXX_FLAGS_DEBUG "-g -static")
set(CMAKE_C_FLAGS_DEBUG "-g -static")

ADD_DEFINITIONS(-D_NO_CEPHES -DOLD_BOOST -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_POSIX_C_SOURCE -D_NOSQLITE)

include_directories (/usr/local/include/boost ${PERCOLATOR_SOURCE_DIR}/src ${PERCOLATOR_SOURCE_DIR}/src/MSToolkit/zLib ${PERCOLATOR_SOURCE_DIR}/src/ArrayLibrary ${PERCOLATOR_SOURCE_DIR}/src/MSToolkit/RAMP ${PERCOLATOR_SOURCE_DIR}/src/MSToolkit/SQLite ${PERCOLATOR_SOURCE_DIR}/src/MSToolkit)

add_subdirectory(ArrayLibrary)
add_subdirectory(MSToolkit)

link_directories(${PERCOLATOR_BINARY_DIR}/src/ArrayLibrary ${PERCOLATOR_BINARY_DIR}/src/MSToolkit)

file(GLOB PERCOLATOR_SOURCES *.cpp)

add_library(perclibrary ${PERCOLATOR_SOURCES})

set(PERC_FILE_NAME percolator)

add_executable(${PERC_FILE_NAME} main.cpp)

target_link_libraries(${PERC_FILE_NAME} perclibrary ArrayLibrary MSToolkit pthread dl)

install(TARGETS ${PERC_FILE_NAME} EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!

add_subdirectory(qvality)
add_subdirectory(elude)
