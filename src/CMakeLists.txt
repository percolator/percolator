set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_C_FLAGS "-O3 ")

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_RELEASE "-O3 ")

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_DEBUG "-g ")

include_directories(${CMAKE_CURRENT_BINARY_DIR})

file(GLOB configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/${myReplFile} @ONLY)
endforeach()


find_package(Xsd REQUIRED)
find_package(XercesC REQUIRED)
find_package(TokyoCabinet REQUIRED)
include_directories(${TOKYOCABINET_INCLUDE_DIR})
include_directories(${XERCESC_INCLUDE_DIR})
include_directories(${XSD_INCLUDE_DIR})

set(xsdfiles  ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.hxx )

# --root-element experiment --root-element frag_spectrum_scan
add_custom_command(
  OUTPUT ${xsdfiles}  
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree --generate-insertion XDR --generate-extraction XDR  --generate-serialization --root-element-all  --namespace-map ${percolator-in-namespace}=percolatorInNs ${CMAKE_CURRENT_SOURCE_DIR}/xml/percolator_in.xsd
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/xml/percolator_in.xsd
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

message("bbbbbbb  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree --generate-serialization --root-element-all --namespace-map http://psidev.info/psi/pi/mzIdentML/1.0=mzIdentML_ns  --namespace-map http://psidev.info/fuge-light/1.0=fuge_light_ns  ${CMAKE_CURRENT_SOURCE_DIR}/mzIdentML1.0.0.xsd ")

set(mzIdentMLxsdfiles ${CMAKE_CURRENT_BINARY_DIR}/mzIdentML1.0.0.cxx ${CMAKE_CURRENT_BINARY_DIR}/mzIdentML1.0.0.hxx )
add_custom_command(
  OUTPUT ${mzIdentMLxsdfiles}  
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree --generate-serialization --root-element-all --namespace-map http://psidev.info/psi/pi/mzIdentML/1.0=mzIdentML_ns  --namespace-map http://psidev.info/fuge-light/1.0=fuge_light_ns  ${CMAKE_CURRENT_SOURCE_DIR}/mzIdentML1.0.0.xsd
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mzIdentML1.0.0.xsd
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

foreach( i ${xsdfiles} ${mzIdentMLxsdfiles} )
  set_source_files_properties(  ${i} PROPERTIES GENERATED true)
  message("dasdfsaas")
endforeach()

ADD_DEFINITIONS(-D_NO_CEPHES -DOLD_BOOST -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_POSIX_C_SOURCE -D_NOSQLITE)

include_directories(. MSToolkit/zLib ArrayLibrary MSToolkit/RAMP MSToolkit)

add_subdirectory(ArrayLibrary)
add_subdirectory(MSToolkit)

file(GLOB PERCOLATOR_SOURCES *.cpp)

add_library(perclibrary ${PERCOLATOR_SOURCES} ${xsdfiles} parser.cxx serializer.cxx)

set(PERC_FILE_NAME percolator)

add_executable(testa test.cpp)
add_executable(${PERC_FILE_NAME} main.cpp)
#add_executable(lightparse lightparse.cpp parser.cxx )

add_executable(convertsequest  ${xsdfiles} ${mzIdentMLxsdfiles} parser.cxx convertsequest.cpp serializer.cxx)


#target_link_libraries(convertsequest perclibrary ${XERCESC_LIBRARIES})
target_link_libraries(convertsequest ${XERCESC_LIBRARIES})


#target_link_libraries(lightparse perclibrary ${XERCESC_LIBRARIES})
target_link_libraries(testa ${XERCESC_LIBRARIES})
target_link_libraries(${PERC_FILE_NAME} perclibrary ArrayLibrary MSToolkit pthread dl ${XERCESC_LIBRARIES} ${TOKYOCABINET_LIBRARIES})

#install(TARGETS lightparse DESTINATION bin)
install(TARGETS convertsequest DESTINATION bin)
install(TARGETS ${PERC_FILE_NAME} EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!

add_subdirectory(qvality)


