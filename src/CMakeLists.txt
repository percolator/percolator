set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_C_FLAGS "-O3 ")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 ")
set(CMAKE_C_FLAGS_RELEASE "-O3 ")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g")


IF(APPLE)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -mmacosx-version-min=10.5")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -mmacosx-version-min=10.5")
  
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mmacosx-version-min=10.5")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mmacosx-version-min=10.5")
  
  set(LINK_FLAGS "-mmacosx-version-min=10.5")
ENDIF(APPLE)

#include_directories (/usr/local/include/boost ${PERCOLATOR_SOURCE_DIR}/src ${PERCOLATOR_SOURCE_DIR}/src/MSToolkit/zLib ${PERCOLATOR_SOURCE_DIR}/src/ArrayLibrary ${PERCOLATOR_SOURCE_DIR}/src/MSToolkit/RAMP ${PERCOLATOR_SOURCE_DIR}/src/MSToolkit/SQLite ${PERCOLATOR_SOURCE_DIR}/src/MSToolkit)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

find_package(Xsd REQUIRED)
find_package(XercesC REQUIRED)

include_directories(${XERCESC_INCLUDE_DIR})
include_directories(${XSD_INCLUDE_DIR})

find_package( Boost 1.36.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(xsdfiles  ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.hxx )
set(percolator_in xml/percolator_in-1-1.xsd)

if(MINGW)
# Installing xdr on mingw is not trivial. Maybe it could be fixed.
  set( xdr_flags  "" )
else()
  set( xdr_flags  --generate-insertion XDR --generate-extraction XDR  )
endif()


add_custom_command(
  OUTPUT ${xsdfiles}   
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree   --generate-serialization ${xdr_flags} --root-element-all  --namespace-map ${percolator-in-namespace}=percolatorInNs ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_in}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_in}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_definitions(-D_NO_CEPHES -DOLD_BOOST -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_POSIX_C_SOURCE -D_NOSQLITE)
add_subdirectory(ArrayLibrary)

include_directories(.)

if(NOT EXCLUDE_CONVERTERS)
  add_subdirectory(converters)
endif()

include_directories(. ArrayLibrary)

# main.cpp

add_library(perclibrary STATIC ${xsdfiles} parser.cxx serializer.cxx  BaseSpline.cpp DescriptionOfCorrect.cpp  MassHandler.cpp  PSMDescription.cpp   ResultHolder.cpp  SetHandler.cpp StdvNormalizer.cpp   svm.cpp   Caller.cpp Enzyme.cpp  Globals.cpp Normalizer.cpp PercolatorCInterface.cpp  SanityCheck.cpp   UniNormalizer.cpp      DataSet.cpp FeatureNames.cpp  LogisticRegression.cpp Option.cpp PosteriorEstimator.cpp    Scores.cpp  SqtSanityCheck.cpp  ssl.cpp EludeModel.cpp )

add_executable(percolator main.cpp)

IF(STATIC)
  set_property(TARGET percolator PROPERTY LINK_SEARCH_END_STATIC ON)
ENDIF()



if(STATIC AND MINGW)
  find_file(XERCESC_STATIC_LIBRARY libxerces-c.a HINTS
     /usr/i686-pc-mingw32/sys-root/mingw/lib/
  )
  message("XERCESC_STATIC_LIBRARY=${XERCESC_STATIC_LIBRARY}")
  set( XERCESC_LIBRARIES ${XERCESC_STATIC_LIBRARY} )
  if(NOT XERCESC_STATIC_LIBRARY)
    message(FATAL_ERROR "could not find the static library libxerces-c.a")
  endif()
  set_target_properties(percolator PROPERTIES LINK_FLAGS "-Wl,-Bstatic -static-libgcc")
endif()


target_link_libraries(percolator perclibrary ArrayLibrary pthread dl ${XERCESC_LIBRARIES})

install(TARGETS percolator EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!

add_subdirectory(qvality)
add_subdirectory(elude)
