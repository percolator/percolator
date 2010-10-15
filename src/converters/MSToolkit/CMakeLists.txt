# Unfortunately there is no standard cmake module for sqlite3
# in http://cmake.org/cmake/help/cmake-2-8-docs.html
# as of 2010-03-05
# Therefore the 
# ${CMAKE_PROJECT_DIR}/cmake/FindSqlite3.cmake
# was added


add_subdirectory(SQLite)
add_subdirectory(zLib)
add_subdirectory(RAMP)

file(GLOB MSTOOLKIT_SOURCES *.cpp)
include_directories(RAMP)
include_directories(Sqlite)
include_directories(zLib)

add_library (MSToolkit ${MSTOOLKIT_SOURCES})

target_link_libraries(MSToolkit zLib SQLite RAMP)



#find_package( Sqlite3 REQUIRED )
#include_directories(${SQLITE3_INCLUDE_DIR})
#
#find_package( ZLIB REQUIRED )
#include_directories(${ZLIB_INCLUDE_DIRS})
#
#add_subdirectory(RAMP)
#
#file(GLOB MSTOOLKIT_SOURCES *.cpp)
#include_directories(RAMP)
#add_library (MSToolkit ${MSTOOLKIT_SOURCES})
#
#target_link_libraries(MSToolkit ${ZLIB_LIBRARIES} ${SQLITE3_LIBRARIES} RAMP )
