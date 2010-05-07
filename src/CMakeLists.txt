
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_C_FLAGS "-O3 ")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 ")
set(CMAKE_C_FLAGS_RELEASE "-O3 ")

set(CMAKE_CXX_FLAGS_DEBUG "-g ")
set(CMAKE_C_FLAGS_DEBUG "-g ")

include_directories(${CMAKE_CURRENT_BINARY_DIR})

file(GLOB configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/${myReplFile} @ONLY)
endforeach()

find_package(Java)
if(NOT JAVA_RUNTIME) 
  message(FATAL_ERROR "Could not find the java runtime")
endif()

find_package(Xsd REQUIRED)
find_package(XercesC REQUIRED)
include_directories(${XERCESC_INCLUDE_DIR})
include_directories(${XSD_INCLUDE_DIR})

#my_find_program(trang)

# todo: search better paths for trang.jar
find_file( TRANG_JAR trang.jar PATHS /afs/pdc.kth.se/home/e/esjolund/Public/trang-20081028/ )
#find_file( TRANG_JAR trang.jar  /tmp/u )


if(NOT TRANG_JAR)
  message(FATAL_ERROR "error: could not find the file trang.jar")
endif()


set(xsdfiles  ${CMAKE_CURRENT_BINARY_DIR}/percolator-xml.cxx )

add_custom_command(
  OUTPUT ${xsdfiles}  ${CMAKE_CURRENT_BINARY_DIR}/percolator-xml.hxx
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree --generate-insertion XDR --generate-extraction XDR  --generate-serialization  --root-element experiment --namespace-map ${percolator-in-namespace}=percolatorInNs ${CMAKE_CURRENT_SOURCE_DIR}/percolator-xml.xsd 
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/percolator-xml.xsd
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
# 


#add_custom_command(
#  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/percolator-xml.xsd 
#  COMMAND ${JAVA_RUNTIME} -jar ${TRANG_JAR} -I rng	 -O xsd  ${CMAKE_CURRENT_BINARY_DIR}/percolator-xml.rng ${CMAKE_CURRENT_BINARY_DIR}/percolator-xml.xsd
#  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/percolator-xml.rng
#  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#)


foreach( i ${xsdfiles} )
  set_source_files_properties(  ${i} PROPERTIES GENERATED true)
endforeach()



ADD_DEFINITIONS(-D_NO_CEPHES -DOLD_BOOST -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_POSIX_C_SOURCE -D_NOSQLITE)

include_directories(. MSToolkit/zLib ArrayLibrary MSToolkit/RAMP MSToolkit)

add_subdirectory(ArrayLibrary)
add_subdirectory(MSToolkit)

#link_directories(${PERCOLATOR_BINARY_DIR}/src/ArrayLibrary ${PERCOLATOR_BINARY_DIR}/src/MSToolkit)

file(GLOB PERCOLATOR_SOURCES *.cpp)

add_library(perclibrary ${PERCOLATOR_SOURCES} ${xsdfiles} parser.cxx)

set(PERC_FILE_NAME percolator)

add_executable(${PERC_FILE_NAME} main.cpp)
add_executable(lightparse lightparse.cpp parser.cxx )
target_link_libraries(lightparse perclibrary ${XERCESC_LIBRARIES})

target_link_libraries(${PERC_FILE_NAME} perclibrary ArrayLibrary MSToolkit pthread dl   ${XERCESC_LIBRARIES})

install(TARGETS ${PERC_FILE_NAME} EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!

add_subdirectory(qvality)




