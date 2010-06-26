set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_C_FLAGS "-O3 ")

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_RELEASE "-O3 ")

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_DEBUG "-g ")


IF(APPLE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -mmacosx-version-min=10.5")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -mmacosx-version-min=10.5")
  
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mmacosx-version-min=10.5")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mmacosx-version-min=10.5")
  
  set(LINK_FLAGS "-mmacosx-version-min=10.5")
ENDIF(APPLE)

#include_directories (/usr/local/include/boost ${PERCOLATOR_SOURCE_DIR}/src ${PERCOLATOR_SOURCE_DIR}/src/MSToolkit/zLib ${PERCOLATOR_SOURCE_DIR}/src/ArrayLibrary ${PERCOLATOR_SOURCE_DIR}/src/MSToolkit/RAMP ${PERCOLATOR_SOURCE_DIR}/src/MSToolkit/SQLite ${PERCOLATOR_SOURCE_DIR}/src/MSToolkit)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

file(GLOB configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cmake )
foreach(i ${configurefiles})
  string(REPLACE ".cmake" "" myReplFile "${i}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${myReplFile}.cmake ${CMAKE_CURRENT_BINARY_DIR}/${myReplFile} @ONLY)
endforeach()

find_package(Xsd REQUIRED)
find_package(XercesC REQUIRED)
include_directories(${XERCESC_INCLUDE_DIR})
include_directories(${XSD_INCLUDE_DIR})

find_package( Boost 1.36.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(xsdfiles  ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.hxx )

# --root-element experiment --root-element frag_spectrum_scan
add_custom_command(
  OUTPUT ${xsdfiles}   
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree --generate-insertion XDR --generate-extraction XDR  --generate-serialization --root-element-all  --namespace-map ${percolator-in-namespace}=percolatorInNs ${CMAKE_CURRENT_SOURCE_DIR}/xml/percolator_in.xsd
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/xml/percolator_in.xsd
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

ADD_DEFINITIONS(-D_NO_CEPHES -DOLD_BOOST -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_POSIX_C_SOURCE -D_NOSQLITE)
add_subdirectory(ArrayLibrary)

include_directories(.)
add_subdirectory(converters)

include_directories(. ArrayLibrary)

# main.cpp

add_library(perclibrary ${xsdfiles} parser.cxx serializer.cxx  BaseSpline.cpp DescriptionOfCorrect.cpp  MassHandler.cpp  PSMDescription.cpp   ResultHolder.cpp  SetHandler.cpp StdvNormalizer.cpp   svm.cpp   Caller.cpp Enzyme.cpp  Globals.cpp Normalizer.cpp PercolatorCInterface.cpp  SanityCheck.cpp   UniNormalizer.cpp      DataSet.cpp FeatureNames.cpp  LogisticRegression.cpp Option.cpp PosteriorEstimator.cpp    Scores.cpp  SqtSanityCheck.cpp  ssl.cpp EludeModel.cpp )

add_executable(percolator main.cpp)

target_link_libraries(percolator perclibrary ArrayLibrary pthread dl ${XERCESC_LIBRARIES})

install(TARGETS percolator EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!

add_subdirectory(qvality)

add_subdirectory(elude)


