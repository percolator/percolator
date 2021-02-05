###############################################################################
# NAME AND VERSION
###############################################################################
cmake_minimum_required (VERSION 2.6)

if(APPLE AND ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER 3.0)
  # Identify AppleClang's COMPILER_ID as AppleClang instead of just Clang
  # (https://gitlab.kitware.com/cmake/cmake/issues/19067)
  cmake_policy(SET CMP0025 NEW)
endif()

project(PERCOLATOR-CONVERTERS)
set(PERCOLATOR_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")
include("${PERCOLATOR_SOURCE_DIR}/CommonCMake.txt")

###############################################################################
# PREPARING TO INSTALL
###############################################################################
IF(APPLE)
  # Fix linking on 10.14+. See https://stackoverflow.com/questions/54068035
  LINK_DIRECTORIES(/usr/local/lib)
ENDIF()

my_set(CMAKE_BUILD_TYPE "Release" "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
my_set(CMAKE_PREFIX_PATH "../" "Default path to packages")
option(XML_SUPPORT "Choose to support xml input (slower compilation)." ON)
if(XML_SUPPORT)
  add_definitions(-DXML_SUPPORT)
endif(XML_SUPPORT)

# PRINT VARIABLES
MESSAGE( STATUS )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS "Building PERCOLATOR-CONVERTERS:" )
MESSAGE( STATUS "change a configuration variable with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "Indicate the type of XML serialization scheme :
		  -DSERIALIZE=\"Boost\" or -DSERIALIZE=\"TokyoCabinet\"
		  or -DSERIALIZE=\"LevelDB\".
		  By setting up the variable SERIALIZE to either Boost
		  or TokyoCabinet or LevelDB you will choose the serialization
		  scheme that will be used to build Converters.
		  Boost serialization option will be used as default if no option is given.")
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "XML_SUPPORT = ${XML_SUPPORT}" )
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
  set(SERDB "Level")
elseif("${SERIALIZE}" STREQUAL "Boost")
  message( STATUS "Using Boost Serialization scheme")
  add_definitions(-D__BOOSTDB__)
  set(BOOSTDB TRUE)
  set(SERDB "Boost")
elseif("${SERIALIZE}" STREQUAL "TokyoCabinet")
  message( STATUS "Using Tokyo Cabinet Serialization scheme")
  add_definitions(-D__TOKYODB__)
  set(TOKYODB TRUE)
  set(SERDB "Tokyo")
else()
  message( STATUS "Using default Serialization scheme : Boost Serialization scheme")
  add_definitions(-D__BOOSTDB__)
  set(BOOSTDB TRUE)
  set(SERDB "Boost")
endif()



MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS )

# Detect machine architecture, on UNIX:
get_arch("${TARGET_ARCH}")

#SET USEFUL VARIABLES TO BE USED IN THE CODE
if(MINGW OR WIN32) # windows
#  set(WRITABLE_DIR "C:\\\\Docume~1\\\\AllUse~1\\\\Documents\\\\converters\\\\")
  set(WRITABLE_DIR "xml\\\\")
  set(MZIDENTML_SCHEMA_LOCATION "xml-mzIdentML-${MZIDENTML_VERSION_MAJOR}-${MZIDENTML_VERSION_MINOR}\\\\")
  set(GAML_TANDEM_SCHEMA_LOCATION "xml-tandem-${TANDEM_VERSION}\\\\")
  set(TANDEM_SCHEMA_LOCATION "xml-tandem-${TANDEM_VERSION}\\\\")
  add_definitions(-D_WIN32)
else() # linux
  set(WRITABLE_DIR "${CMAKE_INSTALL_PREFIX}/share/percolator/")
  set(MZIDENTML_SCHEMA_LOCATION "xml-${MZIDENTML_VERSION_MAJOR}-${MZIDENTML_VERSION_MINOR}/")
  set(GAML_TANDEM_SCHEMA_LOCATION "xml-${TANDEM_VERSION}/")
  set(TANDEM_SCHEMA_LOCATION "xml-${TANDEM_VERSION}/")
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
  message(STATUS "                     Executable: ${XSD_EXECUTABLE}")
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
  if(TokyoCabinet_FOUND)
    message(STATUS "Package Tokyo Cabinet found : ${TokyoCabinet_INCLUDE_DIR}")
  else(TokyoCabinet_FOUND)
    message(FATAL_ERROR "The package Tokyo Cabinet has not been found")
  endif(TokyoCabinet_FOUND)
  include_directories(${TokyoCabinet_INCLUDE_DIR})

  find_package(BZip2 REQUIRED)
  if(BZIP2_FOUND)
    message(STATUS  "Bzip2 found")
  else(BZIP2_FOUND)
    message(FATAL_ERROR "The package Bzip2 has not been found")
  endif()
endif()

if(BOOSTDB)
  find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS serialization filesystem system REQUIRED)
else(BOOSTDB)
  find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS filesystem system REQUIRED)
endif()
include_directories(${Boost_INCLUDE_DIRS})


FIND_PACKAGE ( Threads REQUIRED )

#find_package(Pthreads)
#if(PTHREADS_FOUND)
#  message(STATUS  "Pthreads found")
#else(PTHREADS_FOUND)
#  message(FATAL_ERROR "The package Pthreads has not been not found")
#endif()
#include_directories(${PTHREAD_INCLUDE_DIR})

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

  find_package(ZLIB)
    if(ZLIB_FOUND)
    message(STATUS "Zlib found")
  else(ZLIB_FOUND)
    message(FATAL_ERROR "The package Zlib has not been found")
  endif(ZLIB_FOUND)

  find_package(Iconv REQUIRED)

endif()

if(UNIX)
    set(CMAKE_PREFIX_PATH ${RPCDIR})
    find_package(RPC QUIET)
    set(HAVE_RPC ${RPC_FOUND})
    if(RPC_FOUND)
        message(STATUS "RPC found, will include headers from: ${RPC_INCLUDE_DIR}")
        message(STATUS "and libraries from                  : ${RPC_LIBRARIES}")
        set(LIBRARIES ${LIBRARIES} ${RPC_LIBRARIES})
        include_directories(${RPC_INCLUDE_DIR})
    else(RPC_FOUND)
        message(FATAL_ERROR "RPC support is mandatory.\n"
        "Note that SunRPC has been removed in glibc-2.26 and later, "
        "while being optional in earlier versions. Consider using the "
        "recommended and more modern libtirpc instead.\n"
        "Use -DRPCDIR=DIR to specify the rpc directory tree.\n")
    endif(RPC_FOUND)
endif()


if(WIN32 AND BOOSTDB)
  find_package( ZLIB REQUIRED )
  message(STATUS "ZLIB found: ${ZLIB_INCLUDE_DIRS}")
  if(ZLIB_FOUND)
  else(ZLIB_FOUND)
    message(FATAL_ERROR "The package ZLIB has not been found")
  endif(ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIRS})

  find_package( Sqlite3 REQUIRED )
  if(SQLITE3_FOUND)
    message(STATUS "SQL Lite 3 found: ${SQLITE3_INCLUDE_DIR}")
  else(SQLITE3_FOUND)
    message(FATAL_ERROR "The package SQL Lite 3 has not been found")
  endif(SQLITE3_FOUND)
  include_directories(${SQLITE3_INCLUDE_DIR})
endif()

if(TOKYODB AND APPLE)
  find_package(BZip2 REQUIRED)
endif()



###############################################################################
# RUN CODESYNTHESIS
###############################################################################

if(BOOSTDB AND NOT TOKYODB AND NOT LEVELDB)
  set( xdr_flags --generate-ostream --hxx-prologue-file ${CMAKE_CURRENT_SOURCE_DIR}/library-prologue.hxx --generate-insertion boost::archive::binary_oarchive --generate-extraction boost::archive::binary_iarchive )
else()
  set( xdr_flags --generate-insertion XDR --generate-extraction XDR  )
endif()


message(STATUS "Generating classes from xsd-schema using CodeSynthesis")
# RUN CODESYNTHESIS TO GENERATE xsdfiles for pin (percolator input)
set(perc_in_xsdfiles ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.hxx )
set(percolator_in xml/percolator_in.xsd)
add_custom_command(
  OUTPUT ${perc_in_xsdfiles}
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree --generate-serialization ${xdr_flags} --root-element-all
  --namespace-map ${percolator-in-namespace}=percolatorInNs ${CMAKE_CURRENT_SOURCE_DIR}/../${percolator_in}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../${percolator_in}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# RUN CODESYNTHESIS TO GENERATE xsdfiles for pout (percolator output)
set(perc_out_xsdfiles ${CMAKE_CURRENT_BINARY_DIR}/percolator_out.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_out.hxx )
set(percolator_out xml/percolator_out.xsd)
add_custom_command(
  OUTPUT ${perc_out_xsdfiles}
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree --generate-serialization ${xdr_flags} --root-element-all
  --namespace-map ${percolator-out-namespace}=percolatorOutNs ${CMAKE_CURRENT_SOURCE_DIR}/../${percolator_out}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../${percolator_out}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

ADD_CUSTOM_TARGET(generate_perc_xsdfiles DEPENDS ${perc_in_xsdfiles} ${perc_out_xsdfiles})


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
message(STATUS "Generating ${gaml_tandemxsdfiles}")
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

foreach( i ${perc_in_xsdfiles} ${perc_out_xsdfiles} ${mzIdentMLxsdfiles} ${gaml_tandemxsdfiles} ${tandemxsdfiles})
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
add_library(perclibrary_part STATIC ${perc_in_xsdfiles} ${perc_out_xsdfiles}
	    ../Option.cpp ../Enzyme.cpp ../Globals.cpp ../MassHandler.cpp ../serializer.cxx ../parser.cxx ../Logger.cpp ../MyException.cpp)

# compile converter base files
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

message( STATUS "Using FragSpectrumScanDatabase${SERDB}db.cpp")
add_library(converters STATIC ${mzIdentMLxsdfiles} ${gaml_tandemxsdfiles} ${tandemxsdfiles}
	       Reader.cpp SqtReader.cpp MzidentmlReader.cpp SequestReader.cpp MsgfplusReader.cpp TandemReader.cpp
	       FragSpectrumScanDatabase.cpp Interface.cpp FragSpectrumScanDatabase${SERDB}db.cpp)

ADD_DEPENDENCIES(converters generate_perc_xsdfiles)

#target_link_libraries(converters ${XERCESC_LIBRARIES} ${PTHREADS_LIBRARIES} ${Boost_LIBRARIES})

target_link_libraries(converters ${XERCESC_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${RPC_LIBRARIES})

if(APPLE)
  target_link_libraries(converters ${OpenMP_CXX_LIBRARIES})
endif(APPLE)

if(MINGW)
  if(LEVELDB)
    target_link_libraries(converters ${MINGWLIB} ${XDR_LIBRARIES} ${LDB_LIBRARIES})
  elseif(TOKYODB)
    target_link_libraries(converters ${MINGWLIB} ${XDR_LIBRARIES} ${ZLIB_LIBRARIES} ${MMAN_LIBRARIES} ${PSAPI_LIBRARIES} ${GLOB_LIBRARIES} ${REGEX_LIBRARIES} ${TokyoCabinet_LIBRARIES} ${BZIP2_LIBRARIES})
  elseif(BOOSTDB)
    target_link_libraries(converters ${MINGWLIB})
  endif()
else(MINGW)
  if(LEVELDB)
    target_link_libraries(converters ${LDB_LIBRARIES})
  elseif(TOKYODB)
    target_link_libraries(converters ${ZLIB_LIBRARIES} ${XDR_LIBRARIES} ${TokyoCabinet_LIBRARIES} ${BZIP2_LIBRARIES})
  endif()
endif()

##create the main executables
add_executable(msgf2pin Msgfplus2Pin.cpp )
add_executable(sqt2pin Sqt2Pin.cpp )
add_executable(tandem2pin Tandem2Pin.cpp )
#add_executable(sequest2pin Sequest2Pin.cpp )

target_link_libraries(msgf2pin converters perclibrary_part MSToolkit)
target_link_libraries(sqt2pin converters perclibrary_part MSToolkit)
target_link_libraries(tandem2pin converters perclibrary_part MSToolkit)
#target_link_libraries(sequest2pin converters perclibrary_part MSToolkit )

IF(MINGW)
  set_target_properties(msgf2pin PROPERTIES LINK_FLAGS "-Wl,-Bdynamic -static-libgcc -static-libstdc++")
  set_target_properties(tandem2pin PROPERTIES LINK_FLAGS "-Wl,-Bdynamic -static-libgcc -static-libstdc++")
  set_target_properties(sqt2pin PROPERTIES LINK_FLAGS "-Wl,-Bdynamic -static-libgcc -static-libstdc++")
#  set_target_properties(sequest2pin PROPERTIES LINK_FLAGS "-Wl,-Bdynamic -liconv")
ENDIF()

if(PROFILING AND UNIX)
  set_target_properties(msgf2pin PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
  set_target_properties(sqt2pin PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
  set_target_properties(tandem2pin PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
#  set_target_properties(sequest2pin PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ${PROFILING_COMMON_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 ${PROFILING_COMMON_FLAGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}-O3 ${PROFILING_COMMON_FLAGS}")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 ${PROFILING_COMMON_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -DNDEBUG ${PROFILING_COMMON_FLAGS}")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g -DNDEBUG ${PROFILING_COMMON_FLAGS}")
endif()


if(APPLE)
  set_property(TARGET msgf2pin PROPERTY LINK_SEARCH_START_STATIC ON)
  set_property(TARGET sqt2pin PROPERTY LINK_SEARCH_START_STATIC ON)
  set_property(TARGET tandem2pin PROPERTY LINK_SEARCH_START_STATIC ON)
  set_property(TARGET msgf2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_property(TARGET sqt2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_property(TARGET tandem2pin PROPERTY LINK_SEARCH_END_STATIC ON)
#  set_property(TARGET sequest2pin PROPERTY LINK_SEARCH_END_STATIC ON)
#  set_property(TARGET msgf2pin PROPERTY MACOSX_BUNDLE ON)
#  set_property(TARGET sqt2pin PROPERTY MACOSX_BUNDLE ON)
#  set_property(TARGET tandem2pin PROPERTY MACOSX_BUNDLE ON)
#  set_property(TARGET sequest2pin PROPERTY MACOSX_BUNDLE ON)
endif(APPLE)


#necessary dlls for the windows exe
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



#if(APPLE)
#    install(TARGETS msgf2pin EXPORT msgf2pin  DESTINATION ./bin BUNDLE DESTINATION ../Applications)
#    install(TARGETS sqt2pin EXPORT sqt2pin DESTINATION ./bin BUNDLE DESTINATION ../Applications)
#    install(TARGETS tandem2pin EXPORT tandem2pin  DESTINATION ./bin BUNDLE DESTINATION ../Applications)
#    install(TARGETS sequest2pin EXPORT sequest2pin DESTINATION ./bin BUNDLE DESTINATION ../Applications)
#else(APPLE)
    install(TARGETS msgf2pin EXPORT msgf2pin  DESTINATION bin)
    install(TARGETS sqt2pin EXPORT sqt2pin DESTINATION bin)
    install(TARGETS tandem2pin EXPORT tandem2pin  DESTINATION bin)
#    install(TARGETS sequest2pin EXPORT sequest2pin DESTINATION bin)
#endif(APPLE)

install(FILES ../xml/mzIdentML1.1.0.xsd DESTINATION ${WRITABLE_DIR}/${MZIDENTML_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
install(FILES ../xml/gaml_tandem1.0.xsd DESTINATION ${WRITABLE_DIR}/${GAML_TANDEM_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
install(FILES ../xml/tandem2011.12.01.1.xsd DESTINATION ${WRITABLE_DIR}/${TANDEM_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)


if( MINGW )
    message( STATUS "    Installing system-libraries: MinGW DLLs." )
    set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${MINGW_PATH}/bin/libxerces-c-3-1.dll ${MINGW_PATH}/bin/libbz2-1.dll ${MINGW_PATH}/bin/libsqlite3-0.dll ${MINGW_PATH}/bin/zlib1.dll  )
#    set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${MINGW_PATH}/bin/libgcc_s_sjlj-1.dll ${MINGW_PATH}/bin/libstdc++-6.dll ${MINGW_PATH}/bin/libxerces-c-3-1.dll ${MINGW_PATH}/bin/libboost_system.dll ${MINGW_PATH}/bin/libboost_filesystem.dll ${MINGW_PATH}/bin/libboost_serialization.dll ${MINGW_PATH}/bin/libbz2-1.dll ${MINGW_PATH}/bin/libsqlite3-0.dll  )
    include( InstallRequiredSystemLibraries )
endif( MINGW )

if(MSVC)
  message( STATUS "    Installing system-libraries: Win32 DLLs." )
  # TODO: Some hacking is needed to install the .dll instead of the .lib. CMake plans to include these DLL variables in the future
  string(REPLACE "xerces-c_3.lib" "../bin/xerces-c_3_1.dll;" XERCESC_DLL_LIBRARIES ${XERCESC_LIBRARIES})
  message( STATUS "      Installing Xerces DLL: ${XERCESC_DLL_LIBRARIES}")
  list(REMOVE_ITEM SQLITE3_LIBRARIES debug)
  list(REMOVE_ITEM SQLITE3_LIBRARIES optimized)
  string(REPLACE "sqlite3.lib" "sqlite3.dll;" SQLITE3_DLL_LIBRARIES ${SQLITE3_LIBRARIES})
  message( STATUS "      Installing Sqlite3 DLL: ${SQLITE3_DLL_LIBRARIES}")
  string(REPLACE "zdll.lib" "zlib1.dll" ZLIB_DLL_LIBRARIES ${ZLIB_LIBRARIES}) # win32
  string(REPLACE "lib/zlib.lib" "bin/zlib.dll" ZLIB_DLL_LIBRARIES ${ZLIB_DLL_LIBRARIES}) # win64
  message( STATUS "      Installing Zlib DLL: ${ZLIB_DLL_LIBRARIES}")
  set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${XERCESC_DLL_LIBRARIES} ${ZLIB_DLL_LIBRARIES} ${SQLITE3_DLL_LIBRARIES} )
  include( InstallRequiredSystemLibraries )
endif(MSVC)

###############################################################################
# PACKAGING
###############################################################################

# PACKAGING OPTIONS: GENERAL
set(PACKAGE_NICKNAME "percolator-converters")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Converters include parsers of different database search engines to percolator input format.")
include("${PERCOLATOR_SOURCE_DIR}/CPack.txt")

set(CPACK_RPM_PACKAGE_DEPENDS "libc6, libgcc1")

include(CPack)
