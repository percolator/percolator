###############################################################################
# NAME AND VERSION
###############################################################################
cmake_minimum_required (VERSION 2.6)
project(PERCOLATOR_ELUDE)
set(PERCOLATOR_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")
include("${PERCOLATOR_SOURCE_DIR}/CommonCMake.txt")

###############################################################################
# PREPARING TO INSTALL
###############################################################################

my_set(CMAKE_BUILD_TYPE "Release" "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
my_set(CMAKE_PREFIX_PATH "../" "Default path to packages")

# PRINT VARIBALES
MESSAGE( STATUS )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS "Building ELUDE:" )
MESSAGE( STATUS "change a configuration variable with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "PERCOLATOR_SOURCE_DIR = ${PERCOLATOR_SOURCE_DIR}" )
MESSAGE( STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}" )
MESSAGE( STATUS "GOOGLE_TEST = ${GOOGLE_TEST}" )
MESSAGE( STATUS "GOOGLE_TEST_PATH = ${GOOGLE_TEST_PATH}" )
MESSAGE( STATUS "TARGET_ARCH = ${TARGET_ARCH}" )
MESSAGE( STATUS "PROFILING = ${PROFILING}")
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS )

 
# Detect machine architecture, on UNIX:
get_arch("${TARGET_ARCH}")

# #SET USEFUL VARIABLES TO BE USED IN THE CODE
if(MINGW OR WIN32) # windows
  set(ELUDE_WRITABLE_DIR "C:\\\\Docume~1\\\\AllUse~1\\\\Documents\\\\elude\\\\")
  set(ELUDE_MODELS_PATH "C:\\\\Docume~1\\\\AllUse~1\\\\Documents\\\\elude\\\\models\\\\")
else() # linux
  set(ELUDE_WRITABLE_DIR "${CMAKE_INSTALL_PREFIX}/share/elude/")
  set(ELUDE_MODELS_PATH "${CMAKE_INSTALL_PREFIX}/share/elude/models/")
endif()


# # STORE NEWLY SET VARIABLES IN *.h.cmake FILES
load_config_files("${CMAKE_CURRENT_SOURCE_DIR}/../")
load_config_files("${CMAKE_CURRENT_SOURCE_DIR}")


###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY ELUDE
###############################################################################

find_package(Boost REQUIRED) 
include_directories(${Boost_INCLUDE_DIRS})
#link_directories(${Boost_INCLUDE_DIRS})


###############################################################################
# COMPILING
###############################################################################
# include and link code from percolator
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${PERCOLATOR_SOURCE_DIR}/src)
link_directories(${PERCOLATOR_SOURCE_DIR}/src)

add_library(eludelibrary STATIC RetentionFeatures.cpp DataManager.cpp EludeMain.cpp LibSVRModel.cpp LibsvmWrapper.cpp SVRModel.h RetentionModel.cpp EludeCaller.cpp  
				  LTSRegression.cpp ../svm.cpp ../Normalizer.cpp ../UniNormalizer.cpp ../StdvNormalizer.cpp 
				  ../Option.cpp ../Enzyme.cpp ../PSMDescription.cpp ../Globals.cpp ../Logger.cpp ../MyException.cpp)

add_executable(elude EludeCaller.cpp)

target_link_libraries(elude eludelibrary ${Boost_LIBRARIES})

IF(MINGW OR WIN32)
  set_property(TARGET elude PROPERTY LINK_SEARCH_END_STATIC ON)
  set_target_properties(elude PROPERTIES LINK_FLAGS "-Wl -Bstatic -static-libgcc")
ENDIF()

if(APPLE)
  set_property(TARGET elude PROPERTY LINK_SEARCH_END_STATIC ON)
  set_property(TARGET elude PROPERTY MACOSX_BUNDLE ON) 
endif(APPLE)

if(PROFILING AND UNIX)
  set_target_properties(elude PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}-O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -DNDEBUG -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g -DNDEBUG -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
endif()

###############################################################################
# TESTING
###############################################################################

# Enabling unit level tests (Ctest)
enable_testing()
# Scheduling unit level tests (this will work only under Unix)
if(GOOGLE_TEST)
  add_definitions( -DPATH_TO_DATA="\"\"${PERCOLATOR_SOURCE_DIR}/src/elude_tool/data/elude\"\"" )
  add_definitions( -DPATH_TO_WRITABLE="\"\"/tmp/\"\"" )  
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  include_directories(${PERCOLATOR_SOURCE_DIR}/src/elude_tool/ ${GOOGLE_TEST_PATH}/include/)
  add_executable( runUnitTests data/unit_tests/Unit_tests_Elude_main.cpp )
  target_link_libraries(runUnitTests ${GOOGLE_TEST_PATH}/libgtest.a eludelibrary)
  add_test( runUnitTests runUnitTests )
endif()
# Scheduling system level tests
add_subdirectory(data/system_tests/)


###############################################################################
# CROSS COMPILING AND INSTALLING
###############################################################################

#TODO: make sure the models are included in the windows installer

if(APPLE)
  install(TARGETS elude EXPORT ELUDE DESTINATION ./bin BUNDLE DESTINATION ../Applications) # Important to use relative path here (used by CPack)!
else(APPLE)
  install(TARGETS elude EXPORT ELUDE DESTINATION bin)
endif(APPLE)


if(MINGW)
  add_custom_command( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/elude-${PACKAGE_VERSION}-win32.exe" COMMAND makensis 
		      "${PERCOLATOR_SOURCE_DIR}/cmake/windows_elude/script.nsi" DEPENDS elude) 
  add_custom_target( win32installer DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/elude-${PACKAGE_VERSION}-win32.exe" )
else(MINGW)
  install(DIRECTORY data/elude/models DESTINATION ${ELUDE_MODELS_PATH}.. FILE_PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
endif()

###############################################################################
# PACKAGING
###############################################################################

# PACKAGING OPTIONS: GENERAL
set(PACKAGE_NICKNAME "elude")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Elude software package for prediction of retention times in mass spectorometry experiments.")
include("${PERCOLATOR_SOURCE_DIR}/CPack.txt")
set(CPACK_RPM_PACKAGE_DEPENDS "libc6, libgcc1")

include(CPack)
