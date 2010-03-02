add_subdirectory(zLib)
add_subdirectory(SQLite)
add_subdirectory(RAMP)

file(GLOB MSTOOLKIT_SOURCES *.cpp)
add_library (MSToolkit ${MSTOOLKIT_SOURCES})

target_link_libraries(MSToolkit zLib SQLite RAMP)
