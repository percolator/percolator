#TODO: 
# _CONSIDER PUTTING COMMON ROUTINES IN A SEPARATE FILE (EG: VERSIONING, .H.CMAKE FILE RECURSIVE SUBSTITUTION, USEFUL MACRO, RUN CODESYNTHESIS)...
# UPDATE README FILE AND RELEASE INSTRUCTIONS

###############################################################################
# NAME AND VERSION
###############################################################################
cmake_minimum_required(VERSION 2.8)
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
MESSAGE( STATUS "GOOGLE_TEST_PATH = ${GOOGLE_TEST_PATH}" )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS )

# VERSION AND NAMESPACES - TO BE THE SAME AS PERCOLATOR'S
set(CPACK_PACKAGE_VERSION_MAJOR "unknown")
set(CPACK_PACKAGE_VERSION_MINOR "00")
set(CPACK_PACKAGE_VERSION_PATCH "0")
#TODO: VERSIONING?? schema versions must be the same for converters and percolator
set(PIN_VERSION_MAJOR "1")
set(PIN_VERSION_MINOR "2")
set(POUT_VERSION_MAJOR "1")
set(POUT_VERSION_MINOR "3")

#SET USEFUL VARIABLES TO BE USED IN THE CODE
if(STATIC AND MINGW) # windows
  set(WRITABLE_DIR "C:\\\\Docume~1\\\\AllUse~1\\\\Documents\\\\percolator\\\\")
  set(PIN_SCHEMA_LOCATION "${WRITABLE_DIR}xml-pin-${PIN_VERSION_MAJOR}-${PIN_VERSION_MINOR}\\\\")
  set(TEMP_DIR "%TEMP%\\\\")
  add_definitions(-D_WIN32)
else() # linux
  set(WRITABLE_DIR "${CMAKE_INSTALL_PREFIX}/share/percolator/")
  set(PIN_SCHEMA_LOCATION "${WRITABLE_DIR}xml-pin-${PIN_VERSION_MAJOR}-${PIN_VERSION_MINOR}/")
  set(TEMP_DIR "/tmp/")
endif()
set(PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(percolator-in-namespace http://per-colator.com/percolator_in/${PIN_VERSION_MAJOR}${PIN_VERSION_MINOR})

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

###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY CONVERTERS
###############################################################################

# useful macro
macro(my_find_program progname)
  find_program(${progname}_EXECUTABLE ${progname})
  if(NOT ${progname}_EXECUTABLE)
    message(FATAL_ERROR  "${progname} not found" )
  endif(NOT ${progname}_EXECUTABLE)
endmacro(my_find_program)

set(CMAKE_MODULE_PATH ${PERCOLATOR_SOURCE_DIR}/cmake)

#FIND LIBRARIES
find_package(TokyoCabinet REQUIRED)
include_directories(${TOKYOCABINET_INCLUDE_DIR})
my_find_program(gengetopt)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

find_package(XercesC REQUIRED) # XSD DEPENDS ON XERCESC
include_directories(${XERCESC_INCLUDE_DIR})
find_package(Xsd REQUIRED) # CODESYNTHESIS
include_directories(${XSD_INCLUDE_DIR})
find_package(Boost 1.36.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

###############################################################################
# RUN CODESYNTHESIS AND GENGETOPT
###############################################################################

if(MINGW)
 #Installing xdr on mingw is not trivial. Maybe it could be fixed.
  set( xdr_flags  "" )
else()
  set( xdr_flags  --generate-insertion XDR --generate-extraction XDR  )
endif()

# RUN CODESYNTHESIS on percolator_in.xsd
set(xsdfiles ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.hxx )
set(percolator_in xml/percolator_in.xsd)
add_custom_command(
  OUTPUT ${xsdfiles}
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree   --generate-serialization ${xdr_flags} --root-element-all  --namespace-map ${percolator-in-namespace}=percolatorInNs ${CMAKE_CURRENT_SOURCE_DIR}/../${percolator_in}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../${percolator_in}
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

foreach( i ${xsdfiles} ${mzIdentMLxsdfiles} )
  set_source_files_properties( ${i} PROPERTIES GENERATED true)
endforeach()

# RUN GENGETOPT (define macro to)
macro(gengetoptmacro progname)
add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/${progname}_cmdline.c   ${CMAKE_CURRENT_BINARY_DIR}/${progname}_cmdline.h 
                   COMMAND ${gengetopt_EXECUTABLE} -i ${CMAKE_CURRENT_BINARY_DIR}/${progname}.ggo --file-name=${progname}_cmdline --output-dir=${CMAKE_CURRENT_BINARY_DIR}
                     DEPENDS   ${CMAKE_CURRENT_BINARY_DIR}/${progname}.ggo  
                     VERBATIM)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${progname}_cmdline.c PROPERTIES GENERATED true)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${progname}_cmdline.h PROPERTIES GENERATED true)
endmacro()

###############################################################################
# COMPILING
###############################################################################

#COMPILING (NEEDED SUBSET OF) PERCOLATOR/FIDO CLASSES
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${PERCOLATOR_SOURCE_DIR}/src ${PERCOLATOR_SOURCE_DIR}/src/fido)
link_directories(${PERCOLATOR_SOURCE_DIR}/src)

add_library(perclibrary_part STATIC ${xsdfiles} ../serializer.cxx ../Enzyme.cpp ../Globals.cpp ../MassHandler.cpp ../DataSet.cpp ../FeatureNames.cpp ../parser.cxx ../PosteriorEstimator.cpp ../Scores.cpp ../SetHandler.cpp ../svm.cpp ../BaseSpline.cpp ../LogisticRegression.cpp ../DescriptionOfCorrect.cpp ../PSMDescription.cpp ../ResultHolder.cpp ../Normalizer.cpp ../StdvNormalizer.cpp ../UniNormalizer.cpp ../Option.cpp ../EludeModel.cpp ../fido/PackedVector.cpp ../fido/PackedMatrix.cpp ../fido/Numerical.cpp ../fido/Vector.cpp ../fido/Set.cpp ../fido/)

#COMPILING MSTOOLKIT
add_subdirectory(MSToolkit)
add_library(converters STATIC SqtReader.cpp FragSpectrumScanDatabase.cpp )

#COMPILING MZIDENTML2PIN
gengetoptmacro( mzidentml2pin )
add_executable(mzidentml2pin  ${mzIdentMLxsdfiles} mzidentml2pin.cpp ${CMAKE_CURRENT_BINARY_DIR}/mzidentml2pin_cmdline.c)
target_link_libraries(mzidentml2pin converters perclibrary_part ${XERCESC_LIBRARIES} ${TOKYOCABINET_LIBRARIES})

#COMPILING SQT2PIN
include_directories(. MSToolkit/RAMP MSToolkit )
add_executable(sqt2pin Sqt2Pin.cpp )
target_link_libraries(sqt2pin converters perclibrary_part ${XERCESC_LIBRARIES} ${TOKYOCABINET_LIBRARIES}  MSToolkit RAMP)

IF(STATIC)
  set_property(TARGET mzidentml2pin PROPERTY LINK_SEARCH_END_STATIC ON)
  set_property(TARGET sqt2pin PROPERTY LINK_SEARCH_END_STATIC ON)
ENDIF()

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

install(TARGETS mzidentml2pin DESTINATION ./bin)
install(TARGETS sqt2pin DESTINATION ./bin)

###############################################################################
# PACKAGING
###############################################################################

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
set(TARGET_ARCH "amd64") # Suitable platform architecture for distribution package: 32bit=i386, 64bit=amd64
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
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxerces-c-dev, libboost-dev, gengetopt, libtokyocabinet-dev, zlib1g-dev, percolator (>=2.01)")

# PACKAGING OPTIONS: RPM
set(CPACK_RPM_PACKAGE_LICENSE "Apache 2.0 license")
set(CPACK_RPM_PACKAGE_VENDOR "Stockholms Universitet")
set(CPACK_RPM_PACKAGE_GROUP "Kll-lab")
set(CPACK_RPM_PACKAGE_DEPENDS "libc6, libgcc1")
set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS}, xerces-c-devel, boost-devel, gengetopt, tokyocabinet-devel, zlib-devel, percolator(>=2.01)")

# Use CPack
# the following MUST appear AFTER CPACK variables have been set
include(CPack)

