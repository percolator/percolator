

###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY PERCOLATOR
###############################################################################

# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(${CMAKE_CURRENT_BINARY_DIR})

MESSAGE( STATUS "Checking packages:" )

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

find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS filesystem system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_INCLUDE_DIRS})


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


###############################################################################
# COMPILE INTERNAL LIBRARIES
###############################################################################

#COMPILE FIDO
add_subdirectory(fido)

#INCLUDE PERCOLATOR HEADERS
include_directories(.) 

#INCLUDE FIDO HEADERS FOR PERCLIBRARY 
include_directories(fido)


#############################################################################
# COMPILE PERCOLATOR
###############################################################################

add_library(perclibrary STATIC ${xsdfiles_in} ${xsdfiles_out} parser.cxx serializer.cxx BaseSpline.cpp DescriptionOfCorrect.cpp MassHandler.cpp PSMDescription.cpp 	ResultHolder.cpp 
								  SetHandler.cpp StdvNormalizer.cpp svm.cpp Caller.cpp Enzyme.cpp Globals.cpp Normalizer.cpp PercolatorCInterface.cpp 
								  SanityCheck.cpp UniNormalizer.cpp DataSet.cpp FeatureNames.cpp LogisticRegression.cpp Option.cpp PosteriorEstimator.cpp 
								  ProteinProbEstimator.cpp ProteinFDRestimator.cpp Scores.cpp SqtSanityCheck.cpp ssl.cpp EludeModel.cpp PackedVector.cpp
								  PackedMatrix.cpp Matrix.cpp Logger.cpp MyException.cpp FidoInterface.cpp Protein.cpp )

# COMPILE PERCOLATOR
add_executable(percolator main.cpp)

if(MINGW OR WIN32 OR APPLE)
  set_property(TARGET percolator PROPERTY LINK_SEARCH_END_STATIC TRUE)
  set_property(TARGET percolator MACOSX_BUNDLE ON)
elif(MINGW OR WIN32)
  set_target_properties(percolator PROPERTIES LINK_FLAGS "-Wl,-Bstatic -static-libgcc")
endif()

if(PROFILING AND UNIX)
  set_target_properties(percolator PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}-O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -DNDEBUG -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g -DNDEBUG -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
endif()
 
# SET LIBRARIES FOR PERCOLATOR
target_link_libraries(percolator perclibrary fido pthread ${XERCESC_LIBRARIES} ${Boost_LIBRARIES})
 
# INSTALL PERCOLATOR
IF(NOT MINGW AND NOT WIN32)
  if(APPLE)
     install(TARGETS percolator EXPORT PERCOLATOR DESTINATION ./bin BUNDLE DESTINATION ../Applications)
  else(APPLE)
    install(TARGETS percolator EXPORT PERCOLATOR DESTINATION ./bin) # Important to use relative path here (used by CPack)!
  endif(APPLE)
endif()

###############################################################################
# COMPILE QVALITY
###############################################################################

add_subdirectory(qvality)
