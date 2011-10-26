#TODO: 
# _CONSIDER PUTTING COMMON ROUTINES IN A SEPARATE FILE (EG: VERSIONING, .H.CMAKE FILE RECURSIVE SUBSTITUTION, USEFUL MACRO, RUN CODESYNTHESIS)...
# UPDATE README FILE AND RELEASE INSTRUCTIONS

###############################################################################
# NAME AND VERSION
###############################################################################
cmake_minimum_required (VERSION 2.4)
# We still (try to) support cmake 2.4, but cmake 2.8 policies also work:
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.6)
    cmake_policy(VERSION 2.8)
elseif("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.4)
    cmake_policy(VERSION 2.6)
endif("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.6)

project(PERCOLATOR-CONVERTERS)
set(PERCOLATOR_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")

# PRINT VARIBALES
MESSAGE( STATUS )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS "Building PERCOLATOR-CONVERTERS:" )
MESSAGE( STATUS "change a configuration variable with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "PERCOLATOR_SOURCE_DIR = ${PERCOLATOR_SOURCE_DIR}" )
MESSAGE( STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}" )
MESSAGE( STATUS "TARGET_ARCH = ${TARGET_ARCH}" )
MESSAGE( STATUS "TOOL CHAIN FILE = ${CMAKE_TOOLCHAIN_FILE}")
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS )

# VERSION AND NAMESPACES - TO BE THE SAME AS PERCOLATOR'S
set(CPACK_PACKAGE_VERSION_MAJOR "2")
set(CPACK_PACKAGE_VERSION_MINOR "03")
set(CPACK_PACKAGE_VERSION_PATCH "0")
#TODO: VERSIONING?? schema versions must be the same for converters and percolator
set(PIN_VERSION_MAJOR "1")
set(PIN_VERSION_MINOR "2")
set(POUT_VERSION_MAJOR "1")
set(POUT_VERSION_MINOR "3")

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_C_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
IF(APPLE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -mmacosx-version-min=10.5")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -mmacosx-version-min=10.5")
  
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mmacosx-version-min=10.5")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mmacosx-version-min=10.5")
  
  set(LINK_FLAGS "-mmacosx-version-min=10.5")
ENDIF(APPLE)

IF(WIN32) 
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".dll" ${CMAKE_FIND_LIBRARY_SUFFIXES})
    SET(CMAKE_FIND_LIBRARY_PREFIXES "" "lib") 
    SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/lib/include ${CMAKE_SOURCE_DIR}/lib/include/log4cxx)
    SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${CMAKE_SOURCE_DIR}/lib/lib)
ELSEIF(UNIX)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so" ".la" ${CMAKE_FIND_LIBRARY_SUFFIXES})
    SET(CMAKE_FIND_LIBRARY_PREFIXES "" "lib")
    SET(CMAKE_INCLUDE_PATH /usr/include ${CMAKE_INCLUDE_PATH})
    SET(CMAKE_LIBRARY_PATH /usr/lib /usr/lib64 ${CMAKE_LIBRARY_PATH})
ENDIF(WIN32)

#SET USEFUL VARIABLES TO BE USED IN THE CODE
if(STATIC AND MINGW) # windows
  set(WRITABLE_DIR "C:\\\\Docume~1\\\\AllUse~1\\\\Documents\\\\percolator\\\\")
  set(PIN_SCHEMA_LOCATION "${WRITABLE_DIR}xml-pin-${PIN_VERSION_MAJOR}-${PIN_VERSION_MINOR}\\\\")
  set(POUT_SCHEMA_LOCATION "${WRITABLE_DIR}xml-pout-${POUT_VERSION_MAJOR}-${POUT_VERSION_MINOR}\\\\")
  set(TEMP_DIR "%TEMP%\\\\")
  add_definitions(-D_WIN32)
else() # linux
  set(WRITABLE_DIR "${CMAKE_INSTALL_PREFIX}/share/percolator/")
  set(PIN_SCHEMA_LOCATION "${WRITABLE_DIR}xml-pin-${PIN_VERSION_MAJOR}-${PIN_VERSION_MINOR}/")
  set(POUT_SCHEMA_LOCATION "${WRITABLE_DIR}xml-pout-${POUT_VERSION_MAJOR}-${POUT_VERSION_MINOR}/")
  set(TEMP_DIR "/tmp/")
endif()
set(PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(percolator-in-namespace http://per-colator.com/percolator_in/${PIN_VERSION_MAJOR}${PIN_VERSION_MINOR})
set(percolator-out-namespace http://per-colator.com/percolator_out/${POUT_VERSION_MAJOR}${POUT_VERSION_MINOR})

#TODO: is there a way to do this recursively insteadof doing it twice?
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


# STATIC/DYNAMIC LIBRARY INCLUSION
if(MINGW)
  OPTION(STATIC "build a static binaries. This is mainly intended for the win32-package and the MACOSX-bundle." ON )
else(MINGW)
  OPTION(STATIC "build a static binaries. This is mainly intended for the win32-package and the MACOSX-bundle." OFF )
endif()


###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY CONVERTERS
###############################################################################


set(CMAKE_MODULE_PATH ${PERCOLATOR_SOURCE_DIR}/cmake)

#FIND LIBRARIES
MESSAGE( STATUS "Checking packages:" )

if(MINGW OR WIN32)
    #Look for precompiled Windows libraries 
  find_library(TOKYOCABINET_STATIC_LIBRARY NAMES TokyoCabinet.lib libtokyocabinet.a libtokyocabinet.dll.a HINTS
     /usr/i686-pc-mingw32/sys-root/mingw/lib/
     /usr/i586-mingw32msvc/sys-root/mingw/lib/
     /usr/i586-mingw32msvc/lib/
     /usr/i686-pc-mingw32/lib/
     /mingw/lib/
     ${PROJECT_SOURCE_DIR}/libs/lib
     ${PROJECT_SOURCE_DIR}/libs/dll
  )
  message(STATUS "TOKYOCABINET_STATIC_LIBRARY=${TOKYOCABINET_STATIC_LIBRARY}")
  set( TOKYOCABINET_LIBRARIES ${TOKYOCABINET_STATIC_LIBRARY} )
  if(NOT TOKYOCABINET_LIBRARIES)
    message(FATAL_ERROR "could not find the static library tokyo cabinet")
  endif() 

  find_path(TOKYOCABINET_DIRECTORY_LIBRARY NAMES tcbdb.h HINTS
     /usr/i686-pc-mingw32/sys-root/mingw/include/
     /usr/i586-mingw32msvc/sys-root/mingw/include/
     /usr/i586-mingw32msvc/include/
     /usr/i686-pc-mingw32/include/
     /mingw/include/
     ${PROJECT_SOURCE_DIR}/libs/TokyoCabinet/
     ${PROJECT_SOURCE_DIR}/libs/include/
  )
  message(STATUS "TOKYOCABINET_DIRECTORY_LIBRARY=${TOKYOCABINET_DIRECTORY_LIBRARY}")
  set( TokyoCabinet_INCLUDE_DIR ${TOKYOCABINET_DIRECTORY_LIBRARY} )
  if(NOT TokyoCabinet_INCLUDE_DIR)
    message(FATAL_ERROR "could not find the include dir of tokyo cabinet")
  endif() 
  include_directories(${TokyoCabinet_INCLUDE_DIR})

else(MINGW OR WIN32)
  find_package(TokyoCabinet REQUIRED)
  if(TOKYOCABINET_FOUND)
    message(STATUS "Package Tokyo Cabinet found : ${TokyoCabinet_INCLUDE_DIR}")
  else(TOKYOCABINET_FOUND)
    message(FATAL_ERROR "The package Tokyo Cabinet has not been found, you can download it and install it from : http://fallabs.com/tokyocabinet/ Alternatively you can run this instruction to install it under debian based systems : sudo atp-get install libtokyocabinet1")
  endif(TOKYOCABINET_FOUND)
  include_directories(${TokyoCabinet_INCLUDE_DIR})
endif()

find_package(XercesC REQUIRED) # XSD DEPENDS ON ${XERCESC_FOUND}
if(XERCESC_FOUND)
  message(STATUS "Package Xerces-C++ found: ${XERCESC_INCLUDE_DIR}")
else(XERCESC_FOUND)
  message(FATAL_ERROR "The package Apache Xerces-C++ has not been found, you can download it and install it from : http://xerces.apache.org/xerces-c/download.cgi Alternatively you can run this instruction to install it under debian based systems : sudo atp-get install libxerces-c-dev ")
endif(XERCESC_FOUND)
include_directories(${XERCESC_INCLUDE_DIR})
# if(WIN32)
#   set(CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../libs/xsd-3.3.0-i686/xsd-3.3.0-i686-windows")
# elseif(${TARGET_ARCH} MATCHES "amd64")
#   set(CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../libs/xsd-3.3.0-i686/xsd-3.3.0-x86_64-linux-gnu")  
# else(${TARGET_ARCH} MATCHES "amd64")
#   set(CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../libs/xsd-3.3.0-i686/xsd-3.3.0-i686-linux-gnu")  
# endif()

find_package(Xsd REQUIRED) # CODESYNTHESIS
if(XSD_FOUND)
  message(STATUS "Package CodeSynthesis XSD found: ${XSD_INCLUDE_DIR}")
else(XSD_FOUND)
  message(FATAL_ERROR "The package CodeSynthesis XSD has not been found, you can download it and install it from : http://www.codesynthesis.com/products/xsd/download.xhtml Alternatively you can run this instruction to install it under debian based systems : sudo apt-get install xsdcxx")
endif(XSD_FOUND)  
include_directories(${XSD_INCLUDE_DIR})

find_package(Boost 1.36.0 REQUIRED)
if(Boost_FOUND)
  message(STATUS "Package Boost-${Boost_VERSION} found : ${Boost_INCLUDE_DIR}")
else(Boost_FOUND)
  message(FATAL_ERROR "The package Boost has not been found, you can download it and install it from : http://www.boost.org/users/download/ Alternatively you can run this instruction to install it under debian based systems : sudo apt-get install libboost-dev libboost-doc")
endif(Boost_FOUND)  
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_INCLUDE_DIRS})


if(STATIC AND MINGW)
  #Look for precompiled Windows libraries 
  find_library(PTHREAD_STATIC_LIBRARY NAMES libpthread.a libpthread.dll.a libpthread.la libpthreadGC2.a HINTS
     /usr/i686-pc-mingw32/sys-root/mingw/lib/
     /usr/i586-mingw32msvc/sys-root/mingw/lib/
     /usr/i586-mingw32msvc/lib/
     /usr/i686-pc-mingw32/lib/
     C:/MinGW/lib/
     /mingw/lib   
  )
  message(STATUS "PTHREAD_STATIC_LIBRARY=${PTHREAD_STATIC_LIBRARY}")
  set( PTHREAD_LIBRARIES ${PTHREAD_STATIC_LIBRARY} )
  if(NOT PTHREAD_STATIC_LIBRARY)
    message(FATAL_ERROR "could not find the static library libpthread")
  endif() 

  find_library(PORTABLEXDR_STATIC_LIBRARY NAMES libportablexdr.a libportablexdr.dll.a libportablexdr.la HINTS
     /usr/i686-pc-mingw32/sys-root/mingw/lib/
     /usr/i586-mingw32msvc/sys-root/mingw/lib/
     /usr/i586-mingw32msvc/lib/
     /usr/i686-pc-mingw32/lib/
     C:/MinGW/lib/
     /mingw/lib   
  )
  message(STATUS "PORTABLEXDR_STATIC_LIBRARY=${PORTABLEXDR_STATIC_LIBRARY}")
  set( PORTABLEXDR_LIBRARIES ${PORTABLEXDR_STATIC_LIBRARY} )
  if(NOT PORTABLEXDR_LIBRARIES)
    message(FATAL_ERROR "could not find the static library portableXdr")
  endif() 
endif()
###############################################################################
# RUN CODESYNTHESIS AND GENGETOPT
###############################################################################

if(MINGW)
 #Installing xdr on mingw is not trivial. Maybe it could be fixed.
  set( xdr_flags --generate-insertion XDR --generate-extraction XDR )
else()
  set( xdr_flags  --generate-insertion XDR --generate-extraction XDR  )
endif()

# RUN CODESYNTHESIS TO GENERATE xsdfiles for pin (percolator input)
set(xsdfiles_in ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.hxx )
set(percolator_in xml/percolator_in.xsd)
add_custom_command(
  OUTPUT ${xsdfiles_in}   
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree   --generate-serialization ${xdr_flags} --root-element-all  --namespace-map ${percolator-in-namespace}=percolatorInNs ${CMAKE_CURRENT_SOURCE_DIR}/../${percolator_in}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../${percolator_in}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# RUN CODESYNTHESIS TO GENERATE xsdfiles for pout (percolator output)
set(xsdfiles_out ${CMAKE_CURRENT_BINARY_DIR}/percolator_out.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_out.hxx )
set(percolator_out xml/percolator_out.xsd)
add_custom_command(
  OUTPUT ${xsdfiles_out}   
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree   --generate-serialization ${xdr_flags} --root-element-all  --namespace-map ${percolator-out-namespace}=percolatorOutNs ${CMAKE_CURRENT_SOURCE_DIR}/../${percolator_out}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../${percolator_out}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# RUN CODESYNTHESIS on mzIdentML1.0.0.xsd
set(mzIdentMLxsdfiles ${CMAKE_CURRENT_BINARY_DIR}/mzIdentML1.0.0.cxx ${CMAKE_CURRENT_BINARY_DIR}/mzIdentML1.0.0.hxx )
add_custom_command(
  OUTPUT ${mzIdentMLxsdfiles}  
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree --generate-serialization  --root-element-all --namespace-map http://psidev.info/psi/pi/mzIdentML/1.0=mzIdentML_ns  --namespace-map http://psidev.info/fuge-light/1.0=fuge_light_ns  ${CMAKE_CURRENT_SOURCE_DIR}/mzIdentML1.0.0.xsd
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mzIdentML1.0.0.xsd
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

foreach( i ${xsdfiles_in} ${xsdfiles_out} ${mzIdentMLxsdfiles} )
  set_source_files_properties( ${i} PROPERTIES GENERATED true)
endforeach()


###############################################################################
# COMPILING
###############################################################################
#Cross compilation
# NSIS package options
# Requires NSIS for creating the package
if(MINGW) 
  add_custom_command( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/converters-${PACKAGE_VERSION}-win32.exe" COMMAND makensis "${CMAKE_BINARY_DIR}/cmake/windows_percolator/script.nsi.converter" DEPENDS sqt2pin mzidentml2pin) 
  add_custom_target( win32installer DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/converters-${PACKAGE_VERSION}-win32.exe" )
endif()


ADD_LIBRARY(xerces-c STATIC IMPORTED)
SET_TARGET_PROPERTIES(xerces-c PROPERTIES IMPORTED_LOCATION
   ${XERCESC_LIBRARIES})
ADD_LIBRARY(portablexdr STATIC IMPORTED)
SET_TARGET_PROPERTIES(portablexdr PROPERTIES IMPORTED_LOCATION
   ${PORTABLEXDR_LIBRARIES})

#COMPILING (NEEDED SUBSET OF) PERCOLATOR/FIDO CLASSES
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${PERCOLATOR_SOURCE_DIR}/src ${PERCOLATOR_SOURCE_DIR}/src/fido)
link_directories(${PERCOLATOR_SOURCE_DIR}/src )

add_library(perclibrary_part STATIC ${xsdfiles_in} ${xsdfiles_out} ../serializer.cxx ../Enzyme.cpp ../Globals.cpp ../MassHandler.cpp ../DataSet.cpp ../FeatureNames.cpp ../parser.cxx ../PosteriorEstimator.cpp ../Scores.cpp ../SetHandler.cpp ../svm.cpp ../BaseSpline.cpp ../LogisticRegression.cpp ../DescriptionOfCorrect.cpp ../PSMDescription.cpp ../ResultHolder.cpp ../Normalizer.cpp ../StdvNormalizer.cpp ../UniNormalizer.cpp ../Option.cpp ../EludeModel.cpp ../fido/PackedVector.cpp ../fido/PackedMatrix.cpp ../fido/Numerical.cpp ../fido/Vector.cpp ../fido/Set.cpp ../fido/)

#COMPILING MSTOOLKIT
add_subdirectory(MSToolkit)
add_library(converters STATIC SqtReader.cpp FragSpectrumScanDatabase.cpp )

#COMPILING MZIDENTML2PIN
add_executable(mzidentml2pin  ${mzIdentMLxsdfiles} mzidentml2pin.cpp )
IF(STATIC AND MINGW)
  set_property(TARGET mzidentml2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_target_properties(mzidentml2pin PROPERTIES LINK_FLAGS "-Wl,-Bstatic -static-libgcc")
ENDIF()
target_link_libraries(mzidentml2pin converters perclibrary_part xerces-c portablexdr ${TOKYOCABINET_LIBRARIES} ${PTHREAD_LIBRARIES})

#COMPILING SQT2PIN
include_directories(. MSToolkit/RAMP MSToolkit )
add_executable(sqt2pin Sqt2Pin.cpp )
IF(STATIC AND MINGW)
  set_property(TARGET sqt2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_target_properties(sqt2pin PROPERTIES LINK_FLAGS "-Wl,-Bstatic -static-libgcc")
ENDIF()
target_link_libraries(sqt2pin converters perclibrary_part  MSToolkit RAMP xerces-c portablexdr ${TOKYOCABINET_LIBRARIES} ${PTHREAD_LIBRARIES} )



###############################################################################
# TESTING
###############################################################################

# Enabling unit level tests (Ctest)
enable_testing()
# Scheduling system level tests
add_subdirectory(data/system_tests/converters)
# Scheduling unit level tests
# UNCOMMENT THIS IF UNIT TESTS WILL BE ADDED IN THE FUTURE
#if(GOOGLE_TEST)
#  add_subdirectory(data/unit_tests/converters)
#endif()

###############################################################################
# INSTALLING
###############################################################################

IF(NOT MINGW)
  install(TARGETS mzidentml2pin EXPORT mzidentml2pin  DESTINATION ./bin)
  install(TARGETS sqt2pin EXPORT sqt2pin DESTINATION ./bin)
endif()
###############################################################################
# PACKAGING
###############################################################################
macro(my_set variable default options)
  if(NOT DEFINED ${variable})
    set(${variable} ${default} CACHE STRING ${options} FORCE)
  endif(NOT DEFINED ${variable})
endmacro(my_set)

# PACKAGING OPTIONS: GENERAL
set(PACKAGE_NICKNAME "percolator-converters")
set(CPACK_SET_DESTDIR TRUE) # Enabling absolute paths for CPack (important!)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Software package for the conversion from different search engine result formats to Percolator input format.")
set(CPACK_GENERATOR "DEB;RPM") # Distribution packaging formats ("make package"); RPM packaging can also be included later
set(CPACK_SOURCE_GENERATOR "TGZ") # This file format is used to package source code ("make package_source")
set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_CURRENT_BINARY_DIR};${CMAKE_PROJECT_NAME};ALL;/")
set(CPACK_PACKAGE_VENDOR "Stockholms Universitet")
set(CPACK_PACKAGE_NAME ${PACKAGE_NICKNAME})
#TODO add README file
#set(CPACK_PACKAGE_DESCRIPTION_FILE "${PERCOLATOR_SOURCE_DIR}/src/converters/README")
set(CPACK_RESOURCE_FILE_LICENSE "${PERCOLATOR_SOURCE_DIR}/src/COPYING")
set(PERLUQVAL_VERSION "v${CPACK_PACKAGE_VERSION_MAJOR}_${CPACK_PACKAGE_VERSION_MINOR}")
my_set(TARGET_ARCH "amd64" "set build architecture: options are i386 and amd64")
set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-linux-${TARGET_ARCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME  "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}-src")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PACKAGE_NICKNAME}-${PERLUQVAL_VERSION}")
set(CPACK_STRIP_FILES TRUE)#       NOT SURE IF THIS IS RECURSIVE...

file(GLOB_RECURSE DOT_FILES_BEGIN ".*") # To be ignored by Cpack
file(GLOB_RECURSE TILD_FILES "*~*") # To be ignored by Cpack
# uncomment to print list of files ignored by Cpack
#message("ignoring the following:")
#message("${DOT_FILES_BEGIN}")
#message("${TILD_FILES}")
set(CPACK_SOURCE_IGNORE_FILES "/CVS/;/.svn/;/.swp$/;cscope.*;/.git/;${CMAKE_CURRENT_BINARY_DIR}/;/.bzr/;/.settings/;${DOT_FILES_BEGIN};${TILD_FILES}")

# PACKAGING OPTIONS: DEB
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Lukas KÃ€ll <lukas.kall@cbr.su.se>")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${TARGET_ARCH})
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxerces-c-dev, libboost-dev, gengetopt, libtokyocabinet-dev, zlib1g-dev, libsqlite3-0, percolator (>=2.01)")

# PACKAGING OPTIONS: RPM
set(CPACK_RPM_PACKAGE_LICENSE "Apache 2.0 license")
set(CPACK_RPM_PACKAGE_VENDOR "Stockholms Universitet")
set(CPACK_RPM_PACKAGE_GROUP "Kll-lab")
set(CPACK_RPM_PACKAGE_DEPENDS "libc6, libgcc1")
set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS}, xerces-c-devel, boost-devel, gengetopt, tokyocabinet-devel, zlib-devel, sqlite3-devel, percolator(>=2.01)")



# Use CPack
# the following MUST appear AFTER CPACK variables have been set
include(CPack)

