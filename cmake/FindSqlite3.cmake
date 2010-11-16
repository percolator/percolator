# Unfortunately there is no standard cmake module for sqlite3
# in http://cmake.org/cmake/help/cmake-2-8-docs.html
# as of 2010-03-05
# very spartanic written... but better than nothing

find_path( SQLITE3_INCLUDE_DIR sqlite3.h  )
find_library(SQLITE3_LIBRARIES NAMES sqlite3 )
if( NOT SQLITE3_INCLUDE_DIR )
  message( FATAL_ERROR "no sqlite3.h found")
endif()
if( NOT SQLITE3_LIBRARIES )
  message( FATAL_ERROR "no libsqlite3 found")
endif()

