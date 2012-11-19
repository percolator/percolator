###############################################################################
# NAME AND VERSION
###############################################################################
cmake_minimum_required (VERSION 2.6)
project(PERCOLATOR-CONVERTERS)
set(PERCOLATOR_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")
include("${PERCOLATOR_SOURCE_DIR}/CommonCMake.txt")

###############################################################################
# PREPARING TO INSTALL
###############################################################################

my_set(CMAKE_BUILD_TYPE "Release" "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
my_set(CMAKE_PREFIX_PATH "../" "Default path to packages")

# PRINT VARIABLES
MESSAGE( STATUS )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS "Building PERCOLATOR-CONVERTERS:" )
MESSAGE( STATUS "change a configuration variable with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "Indicate the type of XML serialization scheme : 
		  -DSERIALIZE=\"Boost\" or -DSERIALIZE=\"TokyoCabinet\" or -DSERIALIZE=\"LevelDB\".
		  By setting up the variable SERIALIZE to either Boost or TokyoCabined or LevelDB 
		  you will chose the serialization scheme that will be used to build Converters.
		  Boost serialization option will be used as default if no option is given.")
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "PERCOLATOR_SOURCE_DIR = ${PERCOLATOR_SOURCE_DIR}" )
MESSAGE( STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}" )
MESSAGE( STATUS "TARGET_ARCH = ${TARGET_ARCH}" )
MESSAGE( STATUS "TOOL CHAIN FILE = ${CMAKE_TOOLCHAIN_FILE}")
MESSAGE( STATUS "SERIALIZE = ${SERIALIZE}")
MESSAGE( STATUS "PROFILING = ${PROFILING}")

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

# Detect machine architecture, on UNIX:
get_arch("${TARGET_ARCH}")

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

# STORE NEWLY SET VARIABLES IN *.h.cmake FILES
load_config_files("${CMAKE_CURRENT_SOURCE_DIR}/../")
load_config_files("${CMAKE_CURRENT_SOURCE_DIR}")

###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY CONVERTERS
###############################################################################

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

if(BOOSTDB)
  find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS serialization filesystem system REQUIRED)
else(BOOSTDB)
  find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS filesystem system REQUIRED)
endif()
include_directories(${Boost_INCLUDE_DIRS})
#link_directories(${Boost_INCLUDE_DIRS})

find_package(Pthreads)
if(PTHREADS_FOUND)
  message(STATUS  "Pthreads found")
else(PTHREADS_FOUND)
  message(FATAL_ERROR "The package Pthreads has not been not found")
endif()
include_directories(${PTHREAD_INCLUDE_DIR})

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
# RUN CODESYNTHESIS
###############################################################################

if(BOOSTDB AND NOT TOKYODB AND NOT LEVELDB)
  set( xdr_flags --generate-ostream --hxx-prologue-file ${CMAKE_CURRENT_SOURCE_DIR}/library-prologue.hxx 
       --generate-insertion boost::archive::binary_oarchive --generate-extraction boost::archive::binary_iarchive )
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

# include MSToolkit 
add_subdirectory(MSToolkit)
include_directories(MSToolkit/RAMP MSToolkit )

# include and link code from percolator
include_directories(${PERCOLATOR_SOURCE_DIR}/src)
link_directories(${PERCOLATOR_SOURCE_DIR}/src)
add_library(perclibrary_part STATIC ${xsdfiles_in} ${xsdfiles_out} 
	    ../Option.cpp ../Enzyme.cpp ../Globals.cpp ../serializer.cxx ../parser.cxx ../Logger.cpp ../MyException.cpp)

# compile converter base files
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
if(LEVELDB)
  add_library(converters STATIC ${mzIdentMLxsdfiles} ${gaml_tandemxsdfiles} ${tandemxsdfiles} 
	       Reader.cpp SqtReader.cpp MzidentmlReader.cpp msgfdbReader.cpp tandemReader.cpp 
	       FragSpectrumScanDatabase.cpp FragSpectrumScanDatabaseLeveldb.cpp)
elseif(TOKYODB)
  add_library(converters STATIC ${mzIdentMLxsdfiles} ${gaml_tandemxsdfiles} ${tandemxsdfiles} 
	       Reader.cpp SqtReader.cpp MzidentmlReader.cpp msgfdbReader.cpp tandemReader.cpp 
	       FragSpectrumScanDatabase.cpp FragSpectrumScanDatabaseTokyodb.cpp)
elseif(BOOSTDB)
  add_library(converters STATIC ${mzIdentMLxsdfiles} ${gaml_tandemxsdfiles} ${tandemxsdfiles} 
	       Reader.cpp SqtReader.cpp MzidentmlReader.cpp msgfdbReader.cpp tandemReader.cpp
	       FragSpectrumScanDatabase.cpp FragSpectrumScanDatabaseBoostdb.cpp)
endif()

target_link_libraries(converters ${XERCESC_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES})
if(MINGW)
  if(LEVELDB)
    target_link_libraries(converters ${MINGWLIB} ${XDR_LIBRARIES} ${LDB_LIBRARIES})
  elseif(TOKYODB)
    target_link_libraries(converters ${MINGWLIB} ${XDR_LIBRARIES} ${BZIP2_LIBRARIES} ${ZLIB_LIBRARIES} 
			      ${MMAN_LIBRARIES} ${PSAPI_LIBRARIES} ${GLOB_LIBRARIES} ${REGEX_LIBRARIES} ${TokyoCabinet_LIBRARIES})
  elseif(BOOSTDB)
    target_link_libraries(converters ${MINGWLIB})
  endif()
else(MINGW)
  if(LEVELDB)
    target_link_libraries(converters ${LDB_LIBRARIES})
  elseif(TOKYODB)
    target_link_libraries(converters ${MINGWLIB} ${XDR_LIBRARIES} ${TokyoCabinet_LIBRARIES})
  endif()
endif()

##create the main executables
add_executable(mzidentml2pin mzidentml2pin.cpp )
add_executable(sqt2pin Sqt2Pin.cpp )
add_executable(tandem2pin tandem2pin.cpp )
add_executable(msgfdb2pin msgfdb2Pin.cpp )

target_link_libraries(mzidentml2pin converters perclibrary_part MSToolkit )
target_link_libraries(sqt2pin converters perclibrary_part MSToolkit )
target_link_libraries(tandem2pin converters perclibrary_part MSToolkit )
target_link_libraries(msgfdb2pin converters perclibrary_part MSToolkit )

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


if(PROFILING AND UNIX)
  set_target_properties(mzidentml2pin PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
  set_target_properties(sqt2pin PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
  set_target_properties(tandem2pin PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
  set_target_properties(msgfdb2pin PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}-O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -DNDEBUG -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g -DNDEBUG -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
endif()


if(APPLE)
  set_property(TARGET mzidentml2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_property(TARGET sqt2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_property(TARGET tandem2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_property(TARGET msgfdb2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_property(TARGET mzidentml2pin PROPERTY MACOSX_BUNDLE ON)
  set_property(TARGET sqt2pin PROPERTY MACOSX_BUNDLE ON)
  set_property(TARGET tandem2pin PROPERTY MACOSX_BUNDLE ON)
  set_property(TARGET msgfdb2pin PROPERTY MACOSX_BUNDLE ON) 
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
  
  install(FILES ../xml/mzIdentML1.1.0.xsd DESTINATION ${MZIDENTML_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
  install(FILES ../xml/gaml_tandem1.0.xsd DESTINATION ${GAML_TANDEM_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
  install(FILES ../xml/tandem2011.12.01.1.xsd DESTINATION ${TANDEM_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
elseif(MINGW)
  add_custom_command( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/converters-${PACKAGE_VERSION}-win32.exe" COMMAND makensis 
		      "${PERCOLATOR_SOURCE_DIR}/cmake/windows_converter/script_converters.nsi" DEPENDS sqt2pin mzidentml2pin tandem2pin msgfdb2pin) 
  add_custom_target( win32installer DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/converters-${PACKAGE_VERSION}-win32.exe" )
endif()
  

###############################################################################
# PACKAGING
###############################################################################

# PACKAGING OPTIONS: GENERAL
set(PACKAGE_NICKNAME "percolator-converters")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Converters include parsers of different database search engines to percolator input format.")
include("${PERCOLATOR_SOURCE_DIR}/CPack.txt")

set(CPACK_RPM_PACKAGE_DEPENDS "libc6, libgcc1")

if(TOKYODB)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxerces-c-dev, libboost-dev, libboost-filesystem-dev (>=1.46), libboost-system-dev, 
				      libtokyocabinet-dev, zlib1g-dev, libsqlite3-dev , percolator (>=2.01)")
  set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS}, xerces-c-devel, boost-devel, boost-filesystem (>=1.46), boost-system, tokyocabinet-devel, 
				   zlib-devel, sqlite3-devel, percolator(>=2.01)")
elseif(LEVELDB)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxerces-c-dev, libboost-dev, libboost-filesystem-dev (>=1.46), libboost-system-dev, 
			             leveldb-dev, zlib1g-dev, libsqlite3-dev , percolator (>=2.01)")
  set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS}, xerces-c-devel, boost-devel (>=1.46), boost-filesystem (>=1.46), boost-system, 
				   leveldb-devel, zlib-devel, sqlite3-devel, percolator(>=2.01)")
elseif(BOOSTDB)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxerces-c-dev, libboost-dev,libboost-serialization-dev, libboost-filesystem-dev (>=1.46), libboost-system-dev, 
                                    zlib1g-dev, libsqlite3-dev, percolator (>=2.01)")
  set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS}, xerces-c-devel, boost-devel (>=1.46),boost-serialization, boost-filesystem (>=1.46), boost-system,
                                 zlib-devel, sqlite3-devel, percolator(>=2.01)")
endif()

include(CPack)