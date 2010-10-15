find_package(TokyoCabinet REQUIRED)
include_directories(${TOKYOCABINET_INCLUDE_DIR})

my_find_program(gengetopt)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# reusing from src/CMakeLists.txt. Maybe one can do this recursively instead
file(GLOB configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/${myReplFile} @ONLY)
endforeach()

set(mzIdentMLxsdfiles ${CMAKE_CURRENT_BINARY_DIR}/mzIdentML1.0.0.cxx ${CMAKE_CURRENT_BINARY_DIR}/mzIdentML1.0.0.hxx )
add_custom_command(
  OUTPUT ${mzIdentMLxsdfiles}  
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree --generate-serialization  --root-element-all --namespace-map http://psidev.info/psi/pi/mzIdentML/1.0=mzIdentML_ns  --namespace-map http://psidev.info/fuge-light/1.0=fuge_light_ns  ${CMAKE_CURRENT_SOURCE_DIR}/mzIdentML1.0.0.xsd
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mzIdentML1.0.0.xsd
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

foreach( i ${xsdfiles} ${mzIdentMLxsdfiles} )
  set_source_files_properties( ${i} PROPERTIES GENERATED true)
endforeach()

add_subdirectory(MSToolkit)
add_library(converters STATIC SqtReader.cpp FragSpectrumScanDatabase.cpp )

#add_executable(lightparse lightparse.cpp  ${xsdfiles} ${CMAKE_SOURCE_DIR}/src/parser.cxx )
#target_link_libraries(lightparse ${XERCESC_LIBRARIES})
#install(TARGETS lightparse DESTINATION bin)

macro(gengetoptmacro progname)
add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/${progname}_cmdline.c   ${CMAKE_CURRENT_BINARY_DIR}/${progname}_cmdline.h 
                   COMMAND ${gengetopt_EXECUTABLE} -i ${CMAKE_CURRENT_BINARY_DIR}/${progname}.ggo --file-name=${progname}_cmdline --output-dir=${CMAKE_CURRENT_BINARY_DIR}
                     DEPENDS   ${CMAKE_CURRENT_BINARY_DIR}/${progname}.ggo  
                     VERBATIM)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${progname}_cmdline.c PROPERTIES GENERATED true)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${progname}_cmdline.h PROPERTIES GENERATED true)
endmacro()

gengetoptmacro( mzidentml2pin )
add_executable(mzidentml2pin  ${mzIdentMLxsdfiles} mzidentml2pin.cpp ${CMAKE_CURRENT_BINARY_DIR}/mzidentml2pin_cmdline.c)
target_link_libraries(mzidentml2pin converters perclibrary ArrayLibrary ${XERCESC_LIBRARIES} ${TOKYOCABINET_LIBRARIES})
install(TARGETS mzidentml2pin DESTINATION bin)

# ${mzIdentMLxsdfiles} 
include_directories(. MSToolkit/zLib ArrayLibrary MSToolkit/RAMP MSToolkit)
add_executable(sqt2pin Sqt2Pin.cpp )
target_link_libraries(sqt2pin converters perclibrary ${XERCESC_LIBRARIES} ${TOKYOCABINET_LIBRARIES}  MSToolkit RAMP  ArrayLibrary )
install(TARGETS sqt2pin DESTINATION bin)


IF(STATIC)
  set_property(TARGET mzidentml2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_property(TARGET sqt2pin PROPERTY LINK_SEARCH_END_STATIC ON)
ENDIF()
