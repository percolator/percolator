###############################################################################
# SET COMPILER OPTIONS
###############################################################################

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

###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY PERCOLATOR
###############################################################################

# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(${CMAKE_CURRENT_BINARY_DIR})

MESSAGE( STATUS "Checking packages:" )

find_package(XercesC REQUIRED) # XSD DEPENDS ON ${XERCESC_FOUND}
if(XERCESC_FOUND)
  message(STATUS "Package Xerces-C++ found: ${XERCESC_INCLUDE_DIR}")
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
  
SET (Boost_ADDITIONAL_VERSIONS "1.36" "1.36.0" "1.41" "1.41.0" "1.39" "1.39.0" "1.42.0" "1.42" "1.43.0" 
			       "1.43." "1.44.0" "1.44.0" "1.45.0" "1.45" "1.46.0" "1.46" "1.47.0" "1.47")
find_package(Boost COMPONENTS filesystem system REQUIRED)
if(Boost_FOUND)
  message(STATUS "Package Boost-${Boost_VERSION} found : ${Boost_INCLUDE_DIR}")
else(Boost_FOUND)
  message(FATAL_ERROR "The package Boost has not been found")
endif(Boost_FOUND)  
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_INCLUDE_DIRS})


find_package(Pthreads)
if(PTHREADS_FOUND)
  message(STATUS  "Pthreads found")
else(PTHREADS_FOUND)
  message(FATAL_ERROR "The package Pthreads has not been found")
endif()

###############################################################################
# RUN CODESYNTHESIS
###############################################################################

if(MINGW)
  set( "" )
else()
  set( xdr_flags  --generate-insertion XDR --generate-extraction XDR  )
endif()

# RUN CODESYNTHESIS TO GENERATE xsdfiles for pin (percolator input)
set(xsdfiles_in ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.hxx )
set(percolator_in xml/percolator_in.xsd)
add_custom_command(
  OUTPUT ${xsdfiles_in}   
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree   --generate-serialization  ${xdr_flags} --root-element-all  --namespace-map ${percolator-in-namespace}=percolatorInNs ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_in}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_in}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# RUN CODESYNTHESIS TO GENERATE xsdfiles for pout (percolator output)
set(xsdfiles_out ${CMAKE_CURRENT_BINARY_DIR}/percolator_out.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_out.hxx )
set(percolator_out xml/percolator_out.xsd)
add_custom_command(
  OUTPUT ${xsdfiles_out}   
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree   --generate-serialization ${xdr_flags} --root-element-all  --namespace-map ${percolator-out-namespace}=percolatorOutNs ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_out}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_out}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_definitions(-D_NO_CEPHES -DOLD_BOOST -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_POSIX_C_SOURCE -D_NOSQLITE)


###############################################################################
# COMPILE INTERNAL LIBRARIES
###############################################################################

#COMPILE FIDO
add_subdirectory(fido)


#INCLUDE PERCOLATOR HEADERS
include_directories(.) 

#INCLUDE FIDO HEADERS FOR PERCLIBRARY 
include_directories(. fido)


#############################################################################
# PERCOLATOR
###############################################################################

add_library(perclibrary STATIC ${xsdfiles_in} ${xsdfiles_out} parser.cxx serializer.cxx BaseSpline.cpp DescriptionOfCorrect.cpp MassHandler.cpp PSMDescription.cpp ResultHolder.cpp SetHandler.cpp StdvNormalizer.cpp svm.cpp Caller.cpp Enzyme.cpp Globals.cpp Normalizer.cpp PercolatorCInterface.cpp SanityCheck.cpp UniNormalizer.cpp DataSet.cpp FeatureNames.cpp LogisticRegression.cpp Option.cpp PosteriorEstimator.cpp ProteinProbEstimator.cpp Scores.cpp SqtSanityCheck.cpp ssl.cpp EludeModel.cpp)

# COMPILE PERCOLATOR
add_executable(percolator main.cpp)

IF(MINGW OR WIN32)
  set_property(TARGET percolator PROPERTY LINK_SEARCH_END_STATIC TRUE)
ENDIF()

if(MINGW OR WIN32)
  set_target_properties(percolator PROPERTIES LINK_FLAGS "-Wl,-Bstatic -static-libgcc")
endif()


 
# SET LIBRARIES FOR PERCOLATOR
target_link_libraries(percolator perclibrary fido pthread ${PTHREADS_LIBRARIES} ${XERCESC_LIBRARIES} ${Boost_LIBRARIES})
 
# INSTALL PERCOLATOR
IF(NOT MINGW AND NOT WIN32)
  install(TARGETS percolator EXPORT PERCOLATOR DESTINATION ./bin) # Important to use relative path here (used by CPack)!
endif()

###############################################################################
# COMPILE QUALITY
###############################################################################

add_subdirectory(qvality)
