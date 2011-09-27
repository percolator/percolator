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
SET( BASE "${PROJECT_SOURCE_DIR}/libs")

find_package(XercesC REQUIRED) # XSD DEPENDS ON ${XERCESC_FOUND}
if(XERCESC_FOUND)
  message(STATUS "Package Xerces-C++ found: ${XERCESC_INCLUDE_DIR}")
else(XERCESC_FOUND)
  message(FATAL_ERROR "The package Apache Xerces-C++ has not been found, you can download it and install it from : http://xerces.apache.org/xerces-c/download.cgi Alternatively you can run this instruction to install it under debian based systems : sudo atp-get install libxerces-c-dev ")
endif(XERCESC_FOUND)
include_directories(${XERCESC_INCLUDE_DIR})

find_package(Xsd REQUIRED) # CODESYNTHESIS
if(XSD_FOUND)
  message(STATUS "Package CodeSynthesis XSD found: ${XSD_EXECUTABLE}")
else(XSD_FOUND)
  SET( CS_VERSION "3.3")

  if(MINGW)
    SET( CS_SITE "http://codesynthesis.com/download/xsd/${CS_VERSION}/windows/i686")
    SET( CS_PACK "xsd-${CS_VERSION}.0-i686-windows")
  elseif(${TARGET_ARCH} MATCHES "amd64")
    SET( CS_SITE "http://codesynthesis.com/download/xsd/${CS_VERSION}/linux-gnu/x86_64")
    SET( CS_PACK "xsd-${CS_VERSION}.0-x86_64-linux-gnu")
    else(${TARGET_ARCH} MATCHES "amd64")
      SET( CS_SITE "http://codesynthesis.com/download/xsd/${CS_VERSION}/linux-gnu/i686")
      SET( CS_PACK "xsd-${CS_VERSION}.0-i686-linux-gnu")
  endif(MINGW)

  MESSAGE(STATUS "Downloading XSD...")

  if(MINGW)
    execute_process(COMMAND wget -nc "${CS_SITE}/${CS_PACK}.zip" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
    execute_process(COMMAND unzip "${CS_PACK}.zip" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
    execute_process(COMMAND rm "${CS_PACK}.zip" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
  else(MINGW)
    execute_process(COMMAND wget -nc "${CS_SITE}/${CS_PACK}.tar.bz2" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
    execute_process(COMMAND tar xjf "${CS_PACK}.tar.bz2" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
    execute_process(COMMAND rm "${CS_PACK}.tar.bz2" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
  endif()
  
  
  if(${rv} EQUAL 0) #installed okay
    message(STATUS "XSD has been downloaded succesfully : ${BASE}/${CS_PACK}")
    SET(CMAKE_PREFIX_PATH "${BASE}/${CS_PACK}")
    find_package(Xsd REQUIRED) # CODESYNTHESIS
  else(${rv} EQUAL 0)
    message(FATAL_ERROR "The package CodeSynthesis XSD has not been found, you can download it and install it from : http://www.codesynthesis.com/products/xsd/download.xhtml .Alternatively you can run this instruction to install it under debian based systems : sudo apt-get install xsdcxx ")
  endif(${rv} EQUAL 0)

endif(XSD_FOUND)  
include_directories(${XSD_INCLUDE_DIR})

SET (Boost_ADDITIONAL_VERSIONS "1.36" "1.36.0" "1.41" "1.41.0" "1.39" "1.39.0" "1.42.0" "1.42" "1.43.0" 
			       "1.43." "1.44.0" "1.44.0" "1.45.0" "1.45" "1.46.0" "1.46" "1.47.0" "1.47")
find_package(Boost COMPONENTS filesystem system REQUIRED)
if(Boost_FOUND)
  message(STATUS "Package Boost-${Boost_VERSION} found : ${Boost_INCLUDE_DIR}")
else(Boost_FOUND)
  message(FATAL_ERROR "The package Boost has not been found, you can download it and install it from : http://www.boost.org/users/download/ .Alternatively you can run this instruction to install it under debian based systems : sudo apt-get install libboost-dev libboost-doc")
endif(Boost_FOUND)  
include_directories(${Boost_INCLUDE_DIRS})

if (UNIX)
  find_library (PTHREAD_LIBRARIES pthread PATHS ${CMAKE_SOURCE_DIR}/lib /usr/lib /usr/local/lib)
  if (NOT PTHREAD_LIBRARIES)
    message (SEND_ERROR "PThread library not found")
  else(NOT PTHREAD_LIBRARIES)
    message( STATUS " Package PThread found")
  endif (NOT PTHREAD_LIBRARIES)
endif (UNIX)
mark_as_advanced (PTHREAD_LIBRARIES)

###############################################################################
# RUN CODESYNTHESIS
###############################################################################

if(MINGW)
 #Installing xdr on mingw is not trivial. Maybe it could be fixed.
  set( xdr_flags  "" )
else()
  set( xdr_flags  --generate-insertion XDR --generate-extraction XDR  )
endif()

# RUN CODESYNTHESIS TO GENERATE xsdfiles for pin (percolator input)
set(xsdfiles_in ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.cxx ${CMAKE_CURRENT_BINARY_DIR}/percolator_in.hxx )
set(percolator_in xml/percolator_in.xsd)
add_custom_command(
  OUTPUT ${xsdfiles_in}   
  COMMAND ${XSD_EXECUTABLE} ARGS cxx-tree   --generate-serialization ${xdr_flags} --root-element-all  --namespace-map ${percolator-in-namespace}=percolatorInNs ${CMAKE_CURRENT_SOURCE_DIR}/${percolator_in}
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

IF(STATIC)
  set_property(TARGET percolator PROPERTY LINK_SEARCH_END_STATIC ON)
ENDIF()

if(STATIC AND MINGW)
  #Look for precompiled Windows libraries
#   find_file(XERCESC_STATIC_LIBRARY NAMES libxerces-c.a libxerces-c.la HINTS
#      /usr/i686-pc-mingw32/sys-root/mingw/lib/
#      /usr/i586-mingw32msvc/sys-root/mingw/lib/
#      /usr/i586-mingw32msvc/lib/
#      /usr/i686-pc-mingw32/lib/
#      /mingw/lib/
#   )
#   message(STATUS "XERCESC_STATIC_LIBRARY=${XERCESC_STATIC_LIBRARY}")
#   set( XERCESC_LIBRARIES ${XERCESC_STATIC_LIBRARY} )
#   if(NOT XERCESC_STATIC_LIBRARY)
#     message(FATAL_ERROR "could not find the static library libxerces-c.a")
#   endif()
  
  find_file(PTHREAD_STATIC_LIBRARY NAMES libpthread.a libpthread.la libpthreadGC2.a HINTS
     /usr/i686-pc-mingw32/sys-root/mingw/lib/
     /usr/i586-mingw32msvc/sys-root/mingw/lib/
     /usr/i586-mingw32msvc/lib/
     /usr/i686-pc-mingw32/lib/
     /mingw/lib/
  )
  message(STATUS "PTHREAD_STATIC_LIBRARY=${PTHREAD_STATIC_LIBRARY}")
  set( PTHREAD_LIBRARIES ${PTHREAD_STATIC_LIBRARY} )
  if(NOT PTHREAD_STATIC_LIBRARY)
    message(FATAL_ERROR "could not find the static library libpthread.a")
  endif() 

  set_target_properties(percolator PROPERTIES LINK_FLAGS "-Wl,-Bstatic -static-libgcc")
endif()



# SET LIBRARIES FOR PERCOLATOR
##TODO Probably I have to look for the pre-compiled windows libraries for at least phtread
target_link_libraries(percolator perclibrary fido ${PTHREAD_LIBRARIES} dl ${XERCESC_LIBRARIES} ${Boost_LIBRARIES})

# INSTALL PERCOLATOR
install(TARGETS percolator EXPORT PERCOLATOR DESTINATION ./bin) # Important to use relative path here (used by CPack)!


###############################################################################
# COMPILE QUALITY
###############################################################################

add_subdirectory(qvality)
