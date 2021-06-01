

###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY PERCOLATOR
###############################################################################

# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(${CMAKE_CURRENT_BINARY_DIR})

MESSAGE( STATUS "Checking packages:" )

if(XML_SUPPORT)
  find_package(XercesC REQUIRED) # XSD DEPENDS ON ${XERCESC_FOUND}
  if(XERCESC_FOUND)
    message(STATUS "Found Xerces-C: ${XERCESC_LIBRARIES}")
    message(STATUS "              : ${XERCESC_INCLUDE_DIR}")
    message(STATUS "       Version: ${XERCESC_VERSION}")
  else(XERCESC_FOUND)
    message(FATAL_ERROR "The package Apache Xerces-C++ has not been found")
  endif(XERCESC_FOUND)
  include_directories(${XERCESC_INCLUDE_DIR})

  find_package(Threads REQUIRED) # XERCESC depends on pthread

  if(APPLE OR MINGW)
    find_package(CURL REQUIRED) # XERCESC depends on curl
  endif()

  message(STATUS "Locating CodeSynthesis library with XSDDIR set to $ENV{XSDDIR}")
  find_package(Xsd REQUIRED) # CODESYNTHESIS
  if(XSD_FOUND)
    message(STATUS "Package CodeSynthesis XSD found: ${XSD_INCLUDE_DIR}, ${XSD_EXECUTABLE}")
  else(XSD_FOUND)
    message(FATAL_ERROR "The package CodeSynthesis XSD has not been found")
  endif(XSD_FOUND)
  include_directories(${XSD_INCLUDE_DIR})

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
endif(XML_SUPPORT)

if(APPLE AND OPENMP_FOUND)
  include_directories(${OpenMP_CXX_INCLUDE_DIRS})
endif(APPLE AND OPENMP_FOUND)

# Boost header libraries only for async io header library (asio.hpp)
find_package(Boost ${BOOST_MIN_VERSION} REQUIRED)
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)
if(WIN32)
  add_definitions(-DBOOST_ALL_NO_LIB) # disable autolinking in boost
endif(WIN32)
include_directories(${Boost_INCLUDE_DIRS})

#########################################
# COMPILE BLAS
#########################################
add_subdirectory(blas)
set (BLAS_LIBRARIES ${BLAS_LIBRARIES} blas)

###############################################################################
# RUN CODESYNTHESIS
###############################################################################

if(XML_SUPPORT)
  if(MINGW OR WIN32)
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
endif(XML_SUPPORT)


#############################################################################
# COMPILE PERCOLATOR
###############################################################################

if(XML_SUPPORT)
  add_library(perclibrary STATIC ${xsdfiles_in} ${xsdfiles_out} parser.cxx serializer.cxx BaseSpline.cpp DescriptionOfCorrect.cpp MassHandler.cpp
                  PSMDescription.cpp PSMDescriptionDOC.cpp ResultHolder.cpp
								  XMLInterface.cpp SetHandler.cpp StdvNormalizer.cpp svm.cpp Caller.cpp CrossValidation.cpp Enzyme.cpp Globals.cpp NoNormalizer.cpp Normalizer.cpp
								  SanityCheck.cpp UniNormalizer.cpp DataSet.cpp FeatureNames.cpp LogisticRegression.cpp Option.cpp PosteriorEstimator.cpp
								  ProteinProbEstimator.cpp ProteinFDRestimator.cpp Scores.cpp PseudoRandom.cpp SqtSanityCheck.cpp ssl.cpp EludeModel.cpp PackedVector.cpp
								  PackedMatrix.cpp Matrix.cpp Logger.cpp MyException.cpp FidoInterface.cpp ProteinScoreHolder.cpp PickedProteinInterface.cpp FeatureMemoryPool.cpp GoogleAnalytics.cpp Timer.cpp)
else(XML_SUPPORT)
  add_library(perclibrary STATIC BaseSpline.cpp DescriptionOfCorrect.cpp MassHandler.cpp PSMDescription.cpp PSMDescriptionDOC.cpp ResultHolder.cpp
								  XMLInterface.cpp SetHandler.cpp StdvNormalizer.cpp svm.cpp Caller.cpp CrossValidation.cpp Enzyme.cpp Globals.cpp NoNormalizer.cpp Normalizer.cpp
								  SanityCheck.cpp UniNormalizer.cpp DataSet.cpp FeatureNames.cpp LogisticRegression.cpp Option.cpp PosteriorEstimator.cpp
								  ProteinProbEstimator.cpp ProteinFDRestimator.cpp Scores.cpp PseudoRandom.cpp SqtSanityCheck.cpp ssl.cpp EludeModel.cpp PackedVector.cpp
								  PackedMatrix.cpp Matrix.cpp Logger.cpp MyException.cpp FidoInterface.cpp ProteinScoreHolder.cpp PickedProteinInterface.cpp FeatureMemoryPool.cpp GoogleAnalytics.cpp Timer.cpp)
endif(XML_SUPPORT)


###############################################################################
# COMPILE INTERNAL LIBRARIES
###############################################################################

#COMPILE FIDO
add_subdirectory(fido)

#COMPILE PICKED PROTEIN
add_subdirectory(picked_protein)

#INCLUDE PERCOLATOR HEADERS
include_directories(.)

#INCLUDE FIDO HEADERS FOR PERCLIBRARY
include_directories(fido)

#INCLUDE PICKED PROTEIN HEADERS FOR PERCLIBRARY
include_directories(picked_protein)

################################################################################

# COMPILE PERCOLATOR
if(CRUX)
  ADD_DEFINITIONS(-DCRUX)
  include_directories(${CRUX} ${CRUX}/src ${EXT_BINARY_DIR})
ELSE(CRUX)
  add_executable(percolator main.cpp)

  if(APPLE)
    set_property(TARGET percolator PROPERTY LINK_SEARCH_START_STATIC FALSE)
    set_property(TARGET percolator PROPERTY LINK_SEARCH_END_STATIC TRUE)
  elseif(MINGW)
    set_target_properties(percolator PROPERTIES LINK_FLAGS "-Wl,-Bdynamic -static-libgcc")
  endif()

  if(PROFILING AND UNIX)
    set_target_properties(percolator PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -DNDEBUG ${PROFILING_COMMON_FLAGS}")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g -DNDEBUG ${PROFILING_COMMON_FLAGS}")
  endif()

  # SET LIBRARIES FOR PERCOLATOR
  set(COMMON_LIBRARIES perclibrary ${BLAS_LIBRARIES} fido picked_protein ${XERCESC_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CURL_LIBRARIES} ${RPC_LIBRARIES})
  if(COVERAGE)
    set(COMMON_LIBRARIES -fprofile-arcs ${COMMON_LIBRARIES})
  endif(COVERAGE)
  if(MSVC)
    target_link_libraries(percolator ${COMMON_LIBRARIES})
  else(MSVC)
    if (APPLE)
      target_link_libraries(percolator ${COMMON_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
    else(APPLE)
      target_link_libraries(percolator ${COMMON_LIBRARIES} stdc++)
    endif (APPLE)
  endif()

  # INSTALL PERCOLATOR
  install(TARGETS percolator EXPORT PERCOLATOR DESTINATION bin) # Important to use relative path here (used by CPack)!
ENDIF(CRUX)

###############################################################################
# COMPILE QVALITY
###############################################################################

add_subdirectory(qvality)
