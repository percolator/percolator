set(CMAKE_CXX_FLAGS "-O3 -static")
set(CMAKE_C_FLAGS "-O3 -static")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -static")
set(CMAKE_C_FLAGS_RELEASE "-O3 -static")

set(CMAKE_CXX_FLAGS_DEBUG "-g -static")
set(CMAKE_C_FLAGS_DEBUG "-g -static")

ADD_DEFINITIONS(-DVERSION -D_NO_CEPHES -DOLD_BOOST -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_POSIX_C_SOURCE -D_NOSQLITE)

include_directories (/usr/local/include/boost ${LUMINATOR_SOURCE_DIR}/src ${LUMINATOR_SOURCE_DIR}/src/MSToolkit/zLib ${LUMINATOR_SOURCE_DIR}/src/ArrayLibrary ${LUMINATOR_SOURCE_DIR}/src/MSToolkit/RAMP ${LUMINATOR_SOURCE_DIR}/src/MSToolkit/SQLite ${LUMINATOR_SOURCE_DIR}/src/MSToolkit)

add_subdirectory(ArrayLibrary)
add_subdirectory(MSToolkit)

link_directories(${LUMINATOR_BINARY_DIR}/src/ArrayLibrary ${LUMINATOR_BINARY_DIR}/src/MSToolkit)

file(GLOB LUMINATOR_SOURCES *.cpp)

add_library(perclibrary ${LUMINATOR_SOURCES})

set(PERC_FILE_NAME luminator)

add_executable(${PERC_FILE_NAME} main.cpp)

target_link_libraries(${PERC_FILE_NAME} perclibrary ArrayLibrary MSToolkit pthread dl)

install(TARGETS ${PERC_FILE_NAME} EXPORT ${PROJECT_NAME} DESTINATION bin) # Important to use relative path here (used by CPack)!

add_subdirectory(qvality)
