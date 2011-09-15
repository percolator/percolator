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
# FIND EXTERNAL LIBRARIES REQUIRED BY PERCOLATOR
###############################################################################

# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# (FIND AND) INCLUDE HEADERS OF LIBRARIES

# MACRO(INSTALL_PACKAGE _package _packageName)
#   if(${AUTOMATIC_PACKAGES_INSTALLATION} MATCHES "YES")
#     message(STATUS "${_packageName} not found, proceeding to install it")
#     execute_process(COMMAND id -u OUTPUT_VARIABLE outputRoot)
#     if(${outputRoot} EQUAL 0) # user is Root
#       execute_process(COMMAND apt-get install ${_package} -y OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
#       if(${rv} EQUAL 0) #installed okay
# 	find_package(${_packageName} REQUIRED)
#       else(${rv} EQUAL 0)
# 	message(FATAL_ERROR "There has been a problem installing the package.Error: '${error}")
#       endif(${rv} EQUAL 0)
#     else(${outputRoot} EQUAL 0)
#       message(FATAL_ERROR "You must compile Percolator with sudo rights in order to be able to install the package")
#     endif(${outputRoot} EQUAL 0)
#   endif(${AUTOMATIC_PACKAGES_INSTALLATION} MATCHES "YES")
# ENDMACRO(INSTALL_PACKAGE _package _packageName)

MESSAGE( STATUS "Checking packages:" )
SET( BASE "/tmp")

if(GOOGLE_TEST AND NOT GOOGLE_TEST_PATH) ## DO IT IF ONLY IF GOOGLE TEST IS ACTIVATED AND NOT INSTALLED
  find_package(GoogleTest REQUIRED) #TODO : it does not find the library when it's been installed in  temp
  if(NOT GTEST_FOUND)
    SET( GT_SITE "http://googletest.googlecode.com/files")
    SET( GT_PACK "gtest-1.5.0")

    execute_process(COMMAND wget -nc "${GT_SITE}/${GT_PACK}.tar.gz" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
    execute_process(COMMAND tar xzf "${GT_PACK}.tar.gz" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
    execute_process(COMMAND rm "${GT_PACK}.tar.gz" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)

    if(${rv} EQUAL 0) #installed okay
      message(STATUS "GoogleTest has been downloaded succesfully")
    else(${rv} EQUAL 0)
      message(FATAL_ERROR "There has been a problem downloading the package.Error: '${error}")
    endif(${rv} EQUAL 0)

    message(STATUS "Installing GoogleTest")

    execute_process(COMMAND cmake . WORKING_DIRECTORY "${BASE}/${GT_PACK}" OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
    execute_process(COMMAND make WORKING_DIRECTORY "${BASE}/${GT_PACK}" OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
    execute_process(COMMAND make clean WORKING_DIRECTORY "${BASE}/${GT_PACK}" OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
    #TODO compile it properly and delete it after compiling
    if(${rv} EQUAL 0) #compiled okay
      message(STATUS "GoogleTest has been installed succesfully ${_GTEST_INCLUDE_DIR}")
      find_package(GoogleTest REQUIRED)
      set(GOOGLE_TEST_PATH ${_GTEST_INCLUDE_DIR}) ##set the variable accordingly
    else(${rv} EQUAL 0)
      message(FATAL_ERROR "The package Google Test has not been found, you can download it and install it from :\n
			  http://code.google.com/p/googletest/downloads/list .Alternatively you can run this instruction to install it\n
			  under debian based systems : apt-get install google-gadgets-gtk (you must run it as sudo).")
    endif(${rv} EQUAL 0)

  else(NOT GTEST_FOUND)
    message(STATUS "Package Google Test found: ${_GTEST_INCLUDE_DIR}")
  endif(NOT GTEST_FOUND)

endif(GOOGLE_TEST AND NOT GOOGLE_TEST_PATH)

find_package(XercesC REQUIRED) # XSD DEPENDS ON ${XERCESC_FOUND}
if(XERCESC_FOUND)
  message(STATUS "Package Xerces-C++ found: ${XERCESC_INCLUDE_DIR}")
else(XERCESC_FOUND)
  message(FATAL_ERROR "The package Apache Xerces-C++ has not been found, you can download it and install it from :\n
			http://xerces.apache.org/xerces-c/download.cgi .Alternatively you can run this instruction to install it\n
			under debian based systems : atp-get install libxerces-c-dev (you must run it as sudo).")
endif(XERCESC_FOUND)
include_directories(${XERCESC_INCLUDE_DIR})

find_package(Xsd REQUIRED) # CODESYNTHESIS
if(XSD_FOUND)
  message(STATUS "Package CodeSynthesis XSD found: ${XSD_INCLUDE_DIR}")
else(XSD_FOUND)
  SET( CS_VERSION "3.3")

  if(${TARGET_ARCH} MATCHES "amd64")
    SET( CS_SITE "http://codesynthesis.com/download/xsd/${CS_VERSION}/linux-gnu/x86_64")
    SET( CS_PACK "xsd-${CS_VERSION}.0-x86_64-linux-gnu")
  else(${TARGET_ARCH} MATCHES "amd64")
    SET( CS_SITE "http://codesynthesis.com/download/xsd/${CS_VERSION}/linux-gnu/i686")
    SET( CS_PACK "xsd-${CS_VERSION}.0-i686-linux-gnu")
  endif(${TARGET_ARCH} MATCHES "amd64")

  execute_process(COMMAND wget -nc "${CS_SITE}/${CS_PACK}.tar.bz2" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
  execute_process(COMMAND tar xjf "${CS_PACK}.tar.bz2" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
  execute_process(COMMAND rm "${CS_PACK}.tar.bz2" WORKING_DIRECTORY ${BASE} OUTPUT_VARIABLE output ERROR_VARIABLE error RESULT_VARIABLE rv)
  
  if(${rv} EQUAL 0) #installed okay
    message(STATUS "XSD has been downloaded succesfully")
    find_package(Xsd REQUIRED) # CODESYNTHESIS
  else(${rv} EQUAL 0)
    message(FATAL_ERROR "The package CodeSynthesis XSD has not been found, you can download it and install it from :\n
			http://www.codesynthesis.com/products/xsd/download.xhtml .Alternatively you can run this instruction to install it\n
			under debian based systems : apt-get install xsdcxx (you must run it as sudo).")
  endif(${rv} EQUAL 0)

endif(XSD_FOUND)  
include_directories(${XSD_INCLUDE_DIR})

find_package(Boost 1.36.0 REQUIRED)
if(Boost_FOUND)
  message(STATUS "Package Boost-${Boost_VERSION} found : ${Boost_INCLUDE_DIR}")
else(Boost_FOUND)
  message(FATAL_ERROR "The package Boost has not been found, you can download it and install it from :\n
			http://www.boost.org/users/download/ .Alternatively you can run this instruction to install it\n
			under debian based systems : apt-get install libboost-dev libboost-doc (you must run it as sudo).")
endif(Boost_FOUND)  
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
target_link_libraries(percolator perclibrary fido pthread dl ${XERCESC_LIBRARIES})

# INSTALL PERCOLATOR
install(TARGETS percolator EXPORT PERCOLATOR DESTINATION ./bin) # Important to use relative path here (used by CPack)!


###############################################################################
# COMPILE QUALITY
###############################################################################

add_subdirectory(qvality)
