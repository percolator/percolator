# Unfortunately there is no standard cmake module for sqlite3
# in http://cmake.org/cmake/help/cmake-2-8-docs.html
# as of 2010-03-05
# Therefore the 
# ${CMAKE_PROJECT_DIR}/cmake/FindSqlite3.cmake
# was added

find_package( Sqlite3 REQUIRED )
if(SQLITE3_FOUND)
  message(STATUS "SQL Lite 3 found: ${SQLITE3_INCLUDE_DIR}")
else(SQLITE3_FOUND)
  message(FATAL_ERROR "The package SQL Lite 3 has not been found")
endif(SQLITE3_FOUND)
include_directories(${SQLITE3_INCLUDE_DIR})

find_package( ZLIB REQUIRED )
if(ZLIB_FOUND)
  message(STATUS "Zlib found: ${ZLIB_INCLUDE_DIRS}")
else(ZLIB_FOUND)
  message(FATAL_ERROR "The package Zlib has not been found")
endif(ZLIB_FOUND)
include_directories(${ZLIB_INCLUDE_DIRS})

add_subdirectory(RAMP)

file(GLOB MSTOOLKIT_SOURCES *.cpp)
include_directories(RAMP)
add_library (MSToolkit STATIC  ${MSTOOLKIT_SOURCES})

target_link_libraries(MSToolkit ${ZLIB_LIBRARIES} ${SQLITE3_LIBRARIES} RAMP)
