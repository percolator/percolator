###############################################################################
# NAME AND VERSION
###############################################################################
cmake_minimum_required (VERSION 2.6)
project(PERCOLATOR-CONVERTERS)

set(PERCOLATOR_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")

# PRINT VARIBALES
MESSAGE( STATUS )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS "Indicate the type of XML serialization scheme : 
		  -DSERIALIZE=\"Boost\" or -DSERIALIZE=\"TokyoCabinet\" or -DSERIALIZE=\"LevelDB\".
		  By setting up the variable SERIALIZE to either Boost or TokyoCabined or LevelDB 
		  you will chose the serialization scheme that will be used to build Converters.
		  Boost serialization option will be used as default if no option is given.")
MESSAGE( STATUS "Building PERCOLATOR-CONVERTERS:" )
MESSAGE( STATUS "change a configuration variable with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "PERCOLATOR_SOURCE_DIR = ${PERCOLATOR_SOURCE_DIR}" )
MESSAGE( STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}" )
MESSAGE( STATUS "TARGET_ARCH = ${TARGET_ARCH}" )
MESSAGE( STATUS "TOOL CHAIN FILE = ${CMAKE_TOOLCHAIN_FILE}")
MESSAGE( STATUS "SERIALIZE = ${SERIALIZE}")

if("${SERIALIZE}" STREQUAL "LevelDB")
  message( STATUS "Using LevelDB Serialization scheme")
  add_definitions(-D__LEVELDB__)
  set(LEVELDB TRUE)
elseif("${SERIALIZE}" STREQUAL "Boost")
  message( STATUS "Using Boost Serialization scheme")
  add_definitions(-D__BOOSTDB__)
  set(BOOSTDB TRUE)
elseif("${SERIALIZE}" STREQUAL "TokyoCabinet")
  message( STATUS "Using Tokyo Cabinet Serialization scheme")
  add_definitions(-D__TOKYODB__)
  set(TOKYODB TRUE)
else()
  message( STATUS "Using default Serialization scheme : Boost Serialization scheme")
  add_definitions(-D__BOOSTDB__)
  set(BOOSTDB TRUE)
endif()

MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS )

# VERSION AND NAMESPACES - TO BE THE SAME AS PERCOLATOR'S
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "00")
set(CPACK_PACKAGE_VERSION_PATCH "0")
#TODO: VERSIONING?? schema versions must be the same for converters and percolator
set(PIN_VERSION_MAJOR "1")
set(PIN_VERSION_MINOR "3")
set(POUT_VERSION_MAJOR "1")
set(POUT_VERSION_MINOR "4")
set(MZIDENTML_VERSION_MAJOR "1")
set(MZIDENTML_VERSION_MINOR "1")
set(GAML_TANDEM_VERSION_MAJOR "1")
set(GAML_TANDEM_VERSION_MINOR "0")
set(TANDEM_VERSION "2011.12.01.1")


set(CMAKE_CXX_FLAGS "-O3 -Wno-unused-result -Wno-overflow")
set(CMAKE_C_FLAGS "-O3 -Wno-unused-result -Wno-overflow")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused-result -Wno-overflow -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -Wno-unused-result -Wno-overflow -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wno-unused-result -Wno-overflow")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wno-unused-result -Wno-overflow")
IF(APPLE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -mmacosx-version-min=10.5")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -mmacosx-version-min=10.5")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mmacosx-version-min=10.5")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mmacosx-version-min=10.5")
  set(LINK_FLAGS "-mmacosx-version-min=10.5")
ENDIF(APPLE)


# Detect machine architecture, on UNIX:
if(UNIX AND NOT APPLE AND NOT DEFINED ${TARGET_ARCH})
	execute_process(COMMAND dpkg --print-architecture
		OUTPUT_VARIABLE TARGET_ARCH
                OUTPUT_STRIP_TRAILING_WHITESPACE)
	message(STATUS "Architecture: ${TARGET_ARCH}")
endif(UNIX AND NOT APPLE AND NOT DEFINED ${TARGET_ARCH})
if(APPLE AND NOT DEFINED ${TARGET_ARCH})
	execute_process(COMMAND uname -m
		OUTPUT_VARIABLE TARGET_ARCH
                OUTPUT_STRIP_TRAILING_WHITESPACE)
	message(STATUS "Architecture: ${TARGET_ARCH}")
endif(APPLE AND NOT DEFINED ${TARGET_ARCH})


#SET USEFUL VARIABLES TO BE USED IN THE CODE
if(MINGW OR WIN32) # windows
  set(WRITABLE_DIR "C:\\\\Docume~1\\\\AllUse~1\\\\Documents\\\\converters\\\\")
  set(MZIDENTML_SCHEMA_LOCATION "${WRITABLE_DIR}xml-mzIdentML-${MZIDENTML_VERSION_MAJOR}-${MZIDENTML_VERSION_MINOR}\\\\")
  set(GAML_TANDEM_SCHEMA_LOCATION "${WRITABLE_DIR}xml-gaml_tandem-${GAML_TANDEM_VERSION_MAJOR}-${GAML_TANDEM_VERSION_MINOR}\\\\")
  set(TANDEM_SCHEMA_LOCATION "${WRITABLE_DIR}xml-tandem-${TANDEM_VERSION}\\\\")
  add_definitions(-D_WIN32)
else() # linux
  set(WRITABLE_DIR "${CMAKE_INSTALL_PREFIX}/share/converters/")
  set(MZIDENTML_SCHEMA_LOCATION "${WRITABLE_DIR}xml-${MZIDENTML_VERSION_MAJOR}-${MZIDENTML_VERSION_MINOR}/")
  set(GAML_TANDEM_SCHEMA_LOCATION "${WRITABLE_DIR}xml-${TANDEM_VERSION}/")
  set(TANDEM_SCHEMA_LOCATION "${WRITABLE_DIR}xml-${TANDEM_VERSION}/")
endif()
set(PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(percolator-in-namespace http://per-colator.com/percolator_in/${PIN_VERSION_MAJOR}${PIN_VERSION_MINOR})
set(percolator-out-namespace http://per-colator.com/percolator_out/${POUT_VERSION_MAJOR}${POUT_VERSION_MINOR})
set(mzIdentML-namespace "http://psidev.info/psi/pi/mzIdentML/${MZIDENTML_VERSION_MAJOR}.${MZIDENTML_VERSION_MINOR}")
set(gaml_tandem-namespace "http://www.bioml.com/gaml/")
set(tandem-namespace "http://www.thegpm.org/TANDEM/${TANDEM_VERSION}")

# STORE NEWLY SET VARIABLES IN *.h.cmake FILES
file(GLOB configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/../ ${CMAKE_CURRENT_SOURCE_DIR}/../*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/${myReplFile} @ONLY)
endforeach()

file(GLOB configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/${myReplFile} @ONLY)
endforeach()

file(GLOB configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/windows_converter ${CMAKE_CURRENT_SOURCE_DIR}/cmake/windows_converter/*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/windows_converter/${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/cmake/windows_converter/${myReplFile} @ONLY)
endforeach()


# STATIC/DYNAMIC LIBRARY INCLUSION
if(MINGW OR WIN32 OR APPLE)
  OPTION(STATIC "build a static binaries. This is mainly intended for the win32-package and the MACOSX-bundle." ON )
else(MINGW OR WIN32 OR APPLE)
  OPTION(STATIC "build a static binaries. This is mainly intended for the win32-package and the MACOSX-bundle." OFF )
endif()


###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY CONVERTERS
###############################################################################


set(CMAKE_MODULE_PATH ${PERCOLATOR_SOURCE_DIR}/cmake)

#FIND LIBRARIES
MESSAGE( STATUS "Checking packages.." )

find_package(XercesC REQUIRED) # XSD DEPENDS ON ${XERCESC_FOUND}
if(XERCESC_FOUND)
  message(STATUS "Found Xerces-C: ${XERCESC_LIBRARIES}")
  message(STATUS "              : ${XERCESC_INCLUDE_DIR}")
  message(STATUS "       Version: ${XERCESC_VERSION}")
else(XERCESC_FOUND)
  message(FATAL_ERROR "The package Apache Xerces-C++ has not been found")
endif(XERCESC_FOUND)
include_directories(${XERCESC_INCLUDE_DIR})

find_package(Xsd REQUIRED) # CODESYNTHESIS
if(XSD_FOUND)
  message(STATUS "Package CodeSynthesis XSD found: ${XSD_INCLUDE_DIR}")
else(XSD_FOUND)
  message(FATAL_ERROR "The package CodeSynthesis XSD has not been found")
endif(XSD_FOUND)  
include_directories(${XSD_INCLUDE_DIR})

if(LEVELDB)
  find_package(LDB REQUIRED)
  if(LDB_FOUND)
    message(STATUS "Package LevelDB found : ${LDB_INCLUDE_DIR}")
  else(LDB_FOUND)
    message(FATAL_ERROR "The package LevelDB has not been found")
  endif(LDB_FOUND)
  include_directories(${LDB_INCLUDE_DIR})
elseif(TOKYODB)
  find_package(TokyoCabinet REQUIRED)
  if(TOKYOCABINET_FOUND)
    message(STATUS "Package Tokyo Cabinet found : ${TokyoCabinet_INCLUDE_DIR}")
  else(TOKYOCABINET_FOUND)
    message(FATAL_ERROR "The package Tokyo Cabinet has not been found")
  endif(TOKYOCABINET_FOUND)
  include_directories(${TokyoCabinet_INCLUDE_DIR})
endif()

if(APPLE)
  SET(Boost_USE_STATIC_LIBS   ON)
endif(APPLE)

set(BOOST_MIN_VERSION "1.46.0")
SET (Boost_ADDITIONAL_VERSIONS "1.36" "1.36.0" "1.41" "1.41.0" "1.39" "1.39.0" "1.42.0" "1.42" "1.43.0" 
			       "1.43." "1.44.0" "1.44.0" "1.45.0" "1.45" "1.46.0" "1.46" "1.47.0" "1.47" 
			       "1.48" "1.48.0" "1.49" "1.49.0" "1.50.0" "1.50" "1.51.0" "1.51")
			       
if(BOOSTDB)
  find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS serialization filesystem system REQUIRED)
else(BOOSTDB)
  find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS filesystem system REQUIRED)
endif()
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_INCLUDE_DIRS})

find_package(Pthreads)
if(PTHREADS_FOUND)
  message(STATUS  "Pthreads found")
else(PTHREADS_FOUND)
  message(FATAL_ERROR "The package Pthreads has not been not found")
endif()

if(MINGW AND NOT BOOSTDB)
  find_package(XDR)
  if(XDR_FOUND)
    message(STATUS  "XDR found")
  else(XDR_FOUND)
    message(FATAL_ERROR "The package XDR has not been found")
  endif()
endif()

if(TOKYODB AND MINGW)

  find_package(Glob)
  if(GLOB_FOUND)
    message(STATUS  "Glob found")
  else(GLOB_FOUND)
    message(FATAL_ERROR "The package Glob has not been not found")
  endif()

  find_package(BZip2)
  if(BZIP2_FOUND)
    message(STATUS  "Bzip2 found")
  else(BZIP2_FOUND)
    message(FATAL_ERROR "The package Bzip2 has not been found")
  endif()

  find_package(Regex)
  if(REGEX_FOUND)
    message(STATUS  "Regex found")
  else(REGEX_FOUND)
    message(FATAL_ERROR " The package Regex has not been found")
  endif()

  find_package(ZLIB)
    if(ZLIB_FOUND)
    message(STATUS "Zlib found")
  else(ZLIB_FOUND)
    message(FATAL_ERROR "The package Zlib has not been found")
  endif(ZLIB_FOUND)

  find_package(Mman)
    if(MMAN_FOUND)
    message(STATUS "Mman found")
  else(MMAN_FOUND)
    message(FATAL_ERROR "The package Mman has not been found")
  endif(MMAN_FOUND)

  find_package(PSapi)
    if(PSAPI_FOUND)
    message(STATUS "PSapi found")
  else(PSAPI_FOUND)
    message(FATAL_ERROR "The package PSapi has not been found")
  endif(PSAPI_FOUND)

endif()
###############################################################################
# RUN CODESYNTHESIS AND GENGETOPT
###############################################################################

if(BOOSTDB AND NOT TOKYODB AND NOT LEVELDB)
  set( xdr_flags --generate-ostream --hxx-prologue-file ${CMAKE_CURRENT_SOURCE_DIR}/library-prologue.hxx --generate-insertion boost::archive::binary_oarchive --generate-extraction boost::archive::binary_iarchive )
else()
  set( xdr_flags --generate-insertion XDR --generate-extraction XDR  )
endif()

# RUN CODESYNTHESIS TO GENERATE xsdfiles for pin (percolator input)
set(xsdfiles_in ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.hxx )
set(percolator_in xml/percolator_in.xsd)
add_custom_command(
  OUTPUT ${xsdfiles_in}   
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree   --generate-serialization ${xdr_flags} --root-element-all  
  --namespace-map ${percolator-in-namespace}=percolatorInNs ${CMAKE_CURRENT_SOURCE_DIR}/../${percolator_in}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../${percolator_in}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# RUN CODESYNTHESIS TO GENERATE xsdfiles for pout (percolator output)
set(xsdfiles_out ${CMAKE_CURRENT_BINARY_DIR}/percolator_out.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_out.hxx )
set(percolator_out xml/percolator_out.xsd)
add_custom_command(
  OUTPUT ${xsdfiles_out}   
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree   --generate-serialization ${xdr_flags} --root-element-all  
  --namespace-map ${percolator-out-namespace}=percolatorOutNs ${CMAKE_CURRENT_SOURCE_DIR}/../${percolator_out}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../${percolator_out}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# RUN CODESYNTHESIS on mzIdentML1.1.0.xsd
set(mzIdentMLxsdfiles ${CMAKE_CURRENT_BINARY_DIR}/mzIdentML1.1.0.cxx ${CMAKE_CURRENT_BINARY_DIR}/mzIdentML1.1.0.hxx )
set(mzIdentML_in xml/mzIdentML1.1.0.xsd)
add_custom_command(
  OUTPUT ${mzIdentMLxsdfiles}  
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree --generate-serialization ${xdr_flags}  --root-element-all 
  --namespace-map http://psidev.info/psi/pi/mzIdentML/1.1=mzIdentML_ns  ${CMAKE_CURRENT_SOURCE_DIR}/../${mzIdentML_in}
  #--namespace-map http://psidev.info/fuge-light/1.1=fuge_light_ns
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../${mzIdentML_in}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# RUN CODESYNTHESIS on gaml_tandem.xsd
set(gaml_tandemxsdfiles ${CMAKE_CURRENT_BINARY_DIR}/gaml_tandem1.0.cxx ${CMAKE_CURRENT_BINARY_DIR}/gaml_tandem1.0.hxx )
set(gaml_tandem_in xml/gaml_tandem1.0.xsd)
add_custom_command(
  OUTPUT ${gaml_tandemxsdfiles}  
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree --generate-serialization ${xdr_flags}  --root-element-all 
  --namespace-map http://www.bioml.com/gaml/=gaml_tandem_ns ${CMAKE_CURRENT_SOURCE_DIR}/../${gaml_tandem_in}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../${gaml_tandem_in}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# RUN CODESYNTHESIS on tandem.xsd
set(tandemxsdfiles ${CMAKE_CURRENT_BINARY_DIR}/tandem2011.12.01.1.cxx ${CMAKE_CURRENT_BINARY_DIR}/tandem2011.12.01.1.hxx )
set(tandem_in xml/tandem2011.12.01.1.xsd)
add_custom_command(
  OUTPUT ${tandemxsdfiles}  
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree --generate-serialization ${xdr_flags}  --root-element-all 
  --namespace-map http://www.thegpm.org/TANDEM/2011.12.01.1=tandem_ns
  --namespace-map =tandem_ns
  --namespace-map http://www.bioml.com/gaml/=gaml_tandem_ns ${CMAKE_CURRENT_SOURCE_DIR}/../${tandem_in}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../${tandem_in}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

foreach( i ${xsdfiles_in} ${xsdfiles_out} ${mzIdentMLxsdfiles} ${gaml_tandemxsdfiles} ${tandemxsdfiles})
  set_source_files_properties( ${i} PROPERTIES GENERATED true)
endforeach()


###############################################################################
# COMPILING
###############################################################################
#Cross compilation
# NSIS package options
# Requires NSIS for creating the package
if(MINGW) 
  add_custom_command( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/converters-${PACKAGE_VERSION}-win32.exe" COMMAND makensis "${CMAKE_BINARY_DIR}/cmake/windows_converter/script.nsi" DEPENDS sqt2pin mzidentml2pin tandem2pin msgfdb2pin) 
  add_custom_target( win32installer DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/converters-${PACKAGE_VERSION}-win32.exe" )
endif()

# include and link code from percolator
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${PERCOLATOR_SOURCE_DIR}/src)
link_directories(${PERCOLATOR_SOURCE_DIR}/src)
add_library(perclibrary_part STATIC ${xsdfiles_in} ${xsdfiles_out} ../Option.cpp ../Enzyme.cpp ../Globals.cpp ../serializer.cxx ../parser.cxx ../Logger.cpp)

#COMPILING 
add_subdirectory(MSToolkit)
if(LEVELDB)
  add_library(converters STATIC ${mzIdentMLxsdfiles} ${gaml_tandemxsdfiles} ${tandemxsdfiles} Reader.cpp SqtReader.cpp MzidentmlReader.cpp msgfdbReader.cpp tandemReader.cpp FragSpectrumScanDatabase.cpp FragSpectrumScanDatabaseLeveldb.cpp)
elseif(TOKYODB)
  add_library(converters STATIC ${mzIdentMLxsdfiles} ${gaml_tandemxsdfiles} ${tandemxsdfiles} Reader.cpp SqtReader.cpp MzidentmlReader.cpp msgfdbReader.cpp tandemReader.cpp FragSpectrumScanDatabase.cpp FragSpectrumScanDatabaseTokyodb.cpp)
elseif(BOOSTDB)
  add_library(converters STATIC ${mzIdentMLxsdfiles} ${gaml_tandemxsdfiles} ${tandemxsdfiles} Reader.cpp SqtReader.cpp MzidentmlReader.cpp msgfdbReader.cpp tandemReader.cpp FragSpectrumScanDatabase.cpp FragSpectrumScanDatabaseBoostdb.cpp)
endif()

include_directories(. MSToolkit/RAMP MSToolkit )
add_executable(mzidentml2pin mzidentml2pin.cpp )
add_executable(sqt2pin Sqt2Pin.cpp )
add_executable(tandem2pin tandem2pin.cpp )
add_executable(msgfdb2pin msgfdb2Pin.cpp )
if(MINGW)
  if(LEVELDB)          
    target_link_libraries(mzidentml2pin -L${MINGWLIB} converters perclibrary_part MSToolkit RAMP ${XERCESC_LIBRARIES} ${XDR_LIBRARIES} ${LDB_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES})
    target_link_libraries(sqt2pin -L${MINGWLIB} converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${XDR_LIBRARIES} ${LDB_LIBRARIES}  ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
    target_link_libraries(tandem2pin -L${MINGWLIB} converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${XDR_LIBRARIES} ${LDB_LIBRARIES}  ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
    target_link_libraries(msgfdb2pin -L${MINGWLIB} converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${XDR_LIBRARIES} ${LDB_LIBRARIES}  ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
  elseif(TOKYODB)
    target_link_libraries(mzidentml2pin -L${MINGWLIB} converters perclibrary_part MSToolkit RAMP ${XERCESC_LIBRARIES} ${XDR_LIBRARIES} ${BZIP2_LIBRARIES} ${ZLIB_LIBRARIES} ${MMAN_LIBRARIES} ${PSAPI_LIBRARIES} ${GLOB_LIBRARIES} ${REGEX_LIBRARIES} ${TokyoCabinet_LIBRARIES}  ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES})
    target_link_libraries(sqt2pin -L${MINGWLIB} converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${XDR_LIBRARIES} ${BZIP2_LIBRARIES} ${ZLIB_LIBRARIES} ${MMAN_LIBRARIES} ${PSAPI_LIBRARIES}  ${GLOB_LIBRARIES} ${REGEX_LIBRARIES} ${TokyoCabinet_LIBRARIES}  ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
    target_link_libraries(tandem2pin -L${MINGWLIB} converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${XDR_LIBRARIES} ${BZIP2_LIBRARIES} ${ZLIB_LIBRARIES} ${MMAN_LIBRARIES} ${PSAPI_LIBRARIES}  ${GLOB_LIBRARIES} ${REGEX_LIBRARIES} ${TokyoCabinet_LIBRARIES}  ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
    target_link_libraries(msgfdb2pin -L${MINGWLIB} converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${XDR_LIBRARIES} ${BZIP2_LIBRARIES} ${ZLIB_LIBRARIES} ${MMAN_LIBRARIES} ${PSAPI_LIBRARIES}  ${GLOB_LIBRARIES} ${REGEX_LIBRARIES} ${TokyoCabinet_LIBRARIES}  ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
  elseif(BOOSTDB)
    target_link_libraries(mzidentml2pin -L${MINGWLIB} converters perclibrary_part MSToolkit RAMP ${XERCESC_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES})
    target_link_libraries(sqt2pin -L${MINGWLIB} converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
    target_link_libraries(tandem2pin -L${MINGWLIB} converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
    target_link_libraries(msgfdb2pin -L${MINGWLIB} converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
  endif()
else(MINGW)
  if(LEVELDB)
    target_link_libraries(mzidentml2pin converters perclibrary_part MSToolkit RAMP ${XERCESC_LIBRARIES} ${LDB_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES})
    target_link_libraries(sqt2pin converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${LDB_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
    target_link_libraries(tandem2pin converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${LDB_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
    target_link_libraries(msgfdb2pin converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${LDB_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
  elseif(TOKYODB)
    target_link_libraries(mzidentml2pin converters perclibrary_part MSToolkit RAMP ${XERCESC_LIBRARIES} ${TokyoCabinet_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES})
    target_link_libraries(sqt2pin converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${TokyoCabinet_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
    target_link_libraries(tandem2pin converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${TokyoCabinet_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
    target_link_libraries(msgfdb2pin converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${TokyoCabinet_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
  elseif(BOOSTDB)
    target_link_libraries(mzidentml2pin converters perclibrary_part MSToolkit RAMP ${XERCESC_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES})
    target_link_libraries(sqt2pin converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
    target_link_libraries(tandem2pin converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
    target_link_libraries(msgfdb2pin converters perclibrary_part  MSToolkit RAMP ${XERCESC_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES} )
  endif()
endif()


IF(MINGW OR WIN32)
  set_property(TARGET mzidentml2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_target_properties(mzidentml2pin PROPERTIES LINK_FLAGS "-static-libgcc")
  set_property(TARGET sqt2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_target_properties(sqt2pin PROPERTIES LINK_FLAGS "-static-libgcc")
  set_property(TARGET tandem2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_target_properties(tandem2pin PROPERTIES LINK_FLAGS "-static-libgcc")
  set_property(TARGET msgfdb2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_target_properties(msgfdb2pin PROPERTIES LINK_FLAGS "-static-libgcc")
ENDIF()

if(APPLE)
  set_property(TARGET mzidentml2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_property(TARGET sqt2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_property(TARGET tandem2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_property(TARGET msgfdb2pin PROPERTY LINK_SEARCH_END_STATIC ON)
endif(APPLE)

#necessary dlls for the windows exe
IF(MINGW)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin")

  file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
             FILES_MATCHING PATTERN "boost_serialization*.dll")

  file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
             FILES_MATCHING PATTERN "libxerces-c*.dll")

  file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
             FILES_MATCHING PATTERN "libgcc_s_sjlj*.dll")

  file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
             FILES_MATCHING PATTERN "libstdc++*.dll")

  file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
             FILES_MATCHING PATTERN "libsqlite3*.dll")
             
  file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
             FILES_MATCHING PATTERN "boost_system*.dll")
             
  file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
             FILES_MATCHING PATTERN "boost_filesystem*.dll")

  if(LEVELDB)
    file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/
             FILES_MATCHING PATTERN "boost_interprocess*.dlll")
    file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/
             FILES_MATCHING PATTERN "boost_date_time*.dll")
    file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/
             FILES_MATCHING PATTERN "boost_thread*.dll")
    file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/
             FILES_MATCHING PATTERN "libleveldb*.dll")
    file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/
             FILES_MATCHING PATTERN "libportablexdr*.dll")

             
  elseif(TOKYODB)
    file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/
             FILES_MATCHING PATTERN "zlib*.dll")
    file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/
             FILES_MATCHING PATTERN "libtokyocabinet*.dll")
    file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/
             FILES_MATCHING PATTERN "pthreadGC2*.dll")
    file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/
             FILES_MATCHING PATTERN "libgnurx*.dll")
    file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/
             FILES_MATCHING PATTERN "glob*.dll")
    file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/
             FILES_MATCHING PATTERN "libportablexdr*.dll")
    file(COPY "${MING_PATH}/bin/"
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/
             FILES_MATCHING PATTERN "libbz*.dll")

             
  endif()
endif()
###############################################################################
# TESTING
###############################################################################

# Enabling unit level tests (Ctest)
enable_testing()
# Scheduling system level tests
add_subdirectory(data/system_tests/converters)
# Scheduling unit level tests

###############################################################################
# INSTALLING
###############################################################################

IF(NOT MINGW)
  if(APPLE)
    install(TARGETS mzidentml2pin EXPORT mzidentml2pin  DESTINATION ./bin BUNDLE DESTINATION ../Applications)
    install(TARGETS sqt2pin EXPORT sqt2pin DESTINATION ./bin BUNDLE DESTINATION ../Applications)
    install(TARGETS tandem2pin EXPORT mzidentml2pin  DESTINATION ./bin BUNDLE DESTINATION ../Applications)
    install(TARGETS msgfdb2pin EXPORT sqt2pin DESTINATION ./bin BUNDLE DESTINATION ../Applications)
  else(APPLE)
    install(TARGETS mzidentml2pin EXPORT mzidentml2pin  DESTINATION ./bin)
    install(TARGETS sqt2pin EXPORT sqt2pin DESTINATION ./bin)
    install(TARGETS tandem2pin EXPORT mzidentml2pin  DESTINATION ./bin)
    install(TARGETS msgfdb2pin EXPORT sqt2pin DESTINATION ./bin)
  endif(APPLE)
  
  #install(FILES ../xml/mzIdentML1.0.0.xsd DESTINATION ${MZIDENTML_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
  install(FILES ../xml/mzIdentML1.1.0.xsd DESTINATION ${MZIDENTML_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
  install(FILES ../xml/gaml_tandem1.0.xsd DESTINATION ${GAML_TANDEM_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
  install(FILES ../xml/tandem2011.12.01.1.xsd DESTINATION ${TANDEM_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
  
endif()

###############################################################################
# PACKAGING
###############################################################################
if(NOT MINGW AND NOT WIN32)
macro(my_set variable default options)
  if(NOT DEFINED ${variable})
    set(${variable} ${default} CACHE STRING ${options} FORCE)
  endif(NOT DEFINED ${variable})
endmacro(my_set)

# PACKAGING OPTIONS: GENERAL
set(PACKAGE_NICKNAME "percolator-converters")
set(CPACK_SET_DESTDIR TRUE) # Enabling absolute paths for CPack (important!)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Software package for the conversion from different search engine result formats to Percolator input format.")
set(PERLUQVAL_VERSION "v${CPACK_PACKAGE_VERSION_MAJOR}_${CPACK_PACKAGE_VERSION_MINOR}")
my_set(TARGET_ARCH "noarch" "set build architecture: options are i686, i386 and amd64")

IF(APPLE)
	SET(CPACK_GENERATOR "TGZ;TBZ2;OSXX11")
        set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-mac-${TARGET_ARCH}-${DB}")
        if(${TARGET_ARCH} STREQUAL "x86_64")
	  set(CMAKE_OSX_ARCHITECTURES "ppc;i386;ppc64;x86_64" CACHE STRING "Build architectures for OSX" FORCE)
        else(${TARGET_ARCH} STREQUAL "x86_64")
	  set(CMAKE_OSX_ARCHITECTURES ppc;i386 CACHE STRING "Build architectures for OSX" FORCE)
        endif(${TARGET_ARCH} STREQUAL "x86_64")
        SET(MACOSX_BUNDLE_BUNDLE_NAME Converters)
        SET(MACOSX_BUNDLE_BUNDLE_VERSION ${VERSION})
        SET(MACOSX_BUNDLE_SHORT_VERSION_STRING ${VERSION})
        SET(MACOSX_BUNDLE_LONG_VERSION_STRING "Version ${VERSION}")
elseif(UNIX)
	set(CPACK_GENERATOR "DEB;RPM") # Distribution packaging formats ("make package"); RPM packaging can also be included later
        set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-linux-${TARGET_ARCH}-${DB}")
endif()
set(CPACK_SOURCE_GENERATOR "TGZ") # This file format is used to package source code ("make package_source")
set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_CURRENT_BINARY_DIR};${CMAKE_PROJECT_NAME};ALL;/")
set(CPACK_PACKAGE_VENDOR "Stockholms Universitet")
set(CPACK_PACKAGE_NAME ${PACKAGE_NICKNAME})
#TODO add README file
#set(CPACK_PACKAGE_DESCRIPTION_FILE "${PERCOLATOR_SOURCE_DIR}/src/converters/README")
set(CPACK_RESOURCE_FILE_LICENSE "${PERCOLATOR_SOURCE_DIR}/src/COPYING")
set(CPACK_SOURCE_PACKAGE_FILE_NAME  "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-src")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}")
set(CPACK_STRIP_FILES TRUE)#       NOT SURE IF THIS IS RECURSIVE...

file(GLOB_RECURSE DOT_FILES_BEGIN ".*") # To be ignored by Cpack
file(GLOB_RECURSE TILD_FILES "*~*") # To be ignored by Cpack

set(CPACK_SOURCE_IGNORE_FILES "/CVS/;/.svn/;/.swp$/;cscope.*;/.git/;${CMAKE_CURRENT_BINARY_DIR}/;/.bzr/;/.settings/;${DOT_FILES_BEGIN};${TILD_FILES}")

# PACKAGING OPTIONS: DEB
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Lukas Kall <lukas.kall@cbr.su.se>")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${TARGET_ARCH})
if(TOKYODB)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxerces-c-dev, libboost-dev, libboost-filesystem-dev (>=1.46), libboost-system-dev, 
				      libtokyocabinet-dev, zlib1g-dev, libsqlite3-dev , percolator (>=2.01)")
elseif(LEVELDB)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxerces-c-dev, libboost-dev, libboost-filesystem-dev (>=1.46), libboost-system-dev, 
			             leveldb-dev, zlib1g-dev, libsqlite3-dev , percolator (>=2.01)")
elseif(BOOSTDB)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxerces-c-dev, libboost-dev,libboost-serialization-dev, libboost-filesystem-dev (>=1.46), libboost-system-dev, 
                                    zlib1g-dev, libsqlite3-dev, percolator (>=2.01)")
endif()

# PACKAGING OPTIONS: RPM
set(CPACK_RPM_PACKAGE_LICENSE "Apache 2.0 license")
set(CPACK_RPM_PACKAGE_VENDOR "Stockholms Universitet")
set(CPACK_RPM_PACKAGE_GROUP "Kall-lab")
set(CPACK_RPM_PACKAGE_DEPENDS "libc6, libgcc1")
set(CPACK_RPM_PACKAGE_RELOCATABLE FALSE)
if(TOKYODB)
  set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS}, xerces-c-devel, boost-devel, boost-filesystem (>=1.46), boost-system, tokyocabinet-devel, 
				   zlib-devel, sqlite3-devel, percolator(>=2.01)")
elseif(LEVELDB)
  set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS}, xerces-c-devel, boost-devel (>=1.46), boost-filesystem (>=1.46), boost-system, 
				   leveldb-devel, zlib-devel, sqlite3-devel, percolator(>=2.01)")
elseif(BOOSTDB)
  set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS}, xerces-c-devel, boost-devel (>=1.46),boost-serialization, boost-filesystem (>=1.46), boost-system,
                                 zlib-devel, sqlite3-devel, percolator(>=2.01)")
endif()

include(CPack)
endif()
