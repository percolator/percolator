###############################################################################
# SET COMPILER OPTIONS
###############################################################################

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


###############################################################################
# FIND EXTERNAL LIBRARIES (REQUIRED BY PERCOLATOR, CONVERTERS AND ELUDE)
###############################################################################

# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# (FIND AND) INCLUDE HEADERS OF LIBRARIES
find_package(XercesC REQUIRED) # XSD DEPENDS ON XERCESC
include_directories(${XERCESC_INCLUDE_DIR})
find_package(Xsd REQUIRED) # CODESYNTHESIS
include_directories(${XSD_INCLUDE_DIR})
find_package(Boost 1.36.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})


###############################################################################
# RUN CODESYNTHESIS
###############################################################################

if(MINGW)
 #Installing xdr on mingw is not trivial. Maybe it could be fixed.
  set( xdr_flags  "" )
else()
  set( xdr_flags  --generate-insertion XDR --generate-extraction XDR  )
endif()

set(xsdfiles ${CMAKE_CURRENT_BINARY_DIR}/percolator_in-${SCHEMA_VERSION_MAJOR}-${SCHEMA_VERSION_MINOR}.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_in-${SCHEMA_VERSION_MAJOR}-${SCHEMA_VERSION_MINOR}.hxx )
set(percolator_in xml/percolator_in-${SCHEMA_VERSION_MAJOR}-${SCHEMA_VERSION_MINOR}.xsd)
set(percolator_out xml/percolator_out-${SCHEMA_VERSION_MAJOR}-${SCHEMA_VERSION_MINOR}.xsd)

# RUN CODESYNTHESIS TO GENERATE xsdfiles
add_custom_command(
  OUTPUT ${xsdfiles}   
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree   --generate-serialization ${xdr_flags} --root-element-all  --namespace-map ${percolator-in-namespace}=percolatorInNs ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_in}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_in}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_definitions(-D_NO_CEPHES -DOLD_BOOST -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_POSIX_C_SOURCE -D_NOSQLITE)


###############################################################################
# COMPILE INTERNAL LIBRARIES
###############################################################################

# COMPILE ARRAYLIBRARY
add_subdirectory(ArrayLibrary)


# INCLUDE PERCOLATOR HEADERS FOR CONVERTERS
include_directories(.) 

# COMPILE CONVERTERS
if(NOT EXCLUDE_CONVERTERS)
  add_subdirectory(converters)
endif()



# INCLUDE ARRAYLIBRARY HEADERS FOR PERCLIBRARY 
include_directories(. ArrayLibrary)

# DECLARE LIBRARY PERCLIBRARY
add_library(perclibrary STATIC ${xsdfiles} parser.cxx serializer.cxx BaseSpline.cpp DescriptionOfCorrect.cpp MassHandler.cpp PSMDescription.cpp ResultHolder.cpp SetHandler.cpp StdvNormalizer.cpp svm.cpp Caller.cpp Enzyme.cpp Globals.cpp Normalizer.cpp PercolatorCInterface.cpp SanityCheck.cpp UniNormalizer.cpp DataSet.cpp FeatureNames.cpp LogisticRegression.cpp Option.cpp PosteriorEstimator.cpp Scores.cpp SqtSanityCheck.cpp ssl.cpp EludeModel.cpp)


###############################################################################
# PERCOLATOR
###############################################################################

# COMPILE PERCOLATOR
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



# SET LIBRARIES FOR PERCOLATOR
target_link_libraries(percolator perclibrary ArrayLibrary pthread dl ${XERCESC_LIBRARIES})

# INSTALL PERCOLATOR
install(TARGETS percolator EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!


###############################################################################
# COMPILE QUALITY AND ELUDE
###############################################################################

add_subdirectory(qvality)
if(NOT EXCLUDE_ELUDE)
  add_subdirectory(elude)
endif()
