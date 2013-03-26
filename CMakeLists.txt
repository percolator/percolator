###############################################################################
# NAME AND VERSION
###############################################################################
cmake_minimum_required (VERSION 2.6)
project(PERCOLATOR)
set(PERCOLATOR_SOURCE_DIR ${CMAKE_SOURCE_DIR})
include("CommonCMake.txt")

###############################################################################
# PREPARING TO INSTALL
###############################################################################

my_set(CMAKE_BUILD_TYPE "Debug" "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
my_set(CMAKE_PREFIX_PATH "../" "Default path to packages")

# PRINT VARIBALES TO STDOUT
MESSAGE( STATUS )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS "Building PERCOLATOR:" )
MESSAGE( STATUS "change a configuration variable with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}" )
MESSAGE( STATUS "GOOGLE_TEST = ${GOOGLE_TEST}" )
MESSAGE( STATUS "GOOGLE_TEST_PATH = ${GOOGLE_TEST_PATH}" )
MESSAGE( STATUS "TARGET_ARCH = ${TARGET_ARCH}" )
MESSAGE( STATUS "TOOL CHAIN FILE = ${CMAKE_TOOLCHAIN_FILE}")
MESSAGE( STATUS "PROFILING = ${PROFILING}")
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS )


# Detect machine architecture, on UNIX:
get_arch("${TARGET_ARCH}")

# STORE NEWLY SET VARIABLES IN *.h.cmake FILES
load_config_files("${CMAKE_CURRENT_SOURCE_DIR}/../")
load_config_files("${CMAKE_CURRENT_SOURCE_DIR}")

###############################################################################
# COMPILING
###############################################################################

add_subdirectory(src)

###############################################################################
# TESTING
###############################################################################

# Enabling system level tests (Ctest)
enable_testing()
# Scheduling system level tests
add_subdirectory(data/system_tests/percolator)
# Scheduling unit level tests
if(GOOGLE_TEST)
  add_subdirectory(data/unit_tests/percolator)
endif()

###############################################################################
# INSTALLING
###############################################################################


if( MINGW )
    message( STATUS "    Installing system-libraries: MinGW DLLs." )
    set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${MINGW_PATH}/bin/libgcc_s_sjlj-1.dll ${MINGW_PATH}/bin/libstdc++-6.dll ${MINGW_PATH}/bin/libxerces-c-3-1.dll ${MINGW_PATH}/bin/libboost_system.dll ${MINGW_PATH}/bin/libboost_filesystem.dll)
    include( InstallRequiredSystemLibraries )
endif( MINGW )

install(FILES src/xml/percolator_in.xsd DESTINATION xml/${PIN_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
install(FILES src/xml/percolator_out.xsd DESTINATION xml/${POUT_SCHEMA_LOCATION} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)

# NSIS package options
# Requires NSIS for creating the package 
#if(MINGW)
#  add_custom_command( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/percolator-${PACKAGE_VERSION}-win32.exe" 
#		       COMMAND makensis "${CMAKE_BINARY_DIR}/cmake/windows_percolator/script.nsi" DEPENDS percolator qvality) 
#  add_custom_target( win32installer DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/percolator-${PACKAGE_VERSION}-win32.exe" )
#endif()

#IF(MINGW)
#  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin")
#  file(COPY "${MINGW_PATH}/bin/"
#             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
#             FILES_MATCHING PATTERN "boost_system*.dll")
#  file(COPY "${MINGW_PATH}/bin/"
#             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
#             FILES_MATCHING PATTERN "boost_filesystem*.dll")
#  file(COPY "${MINGW_PATH}/bin/"
#             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
#             FILES_MATCHING PATTERN "libxerces-c*.dll")
#  file(COPY "${MINGW_PATH}/bin/"
#             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
#             FILES_MATCHING PATTERN "libgcc_s_sjlj*.dll")
#  file(COPY "${MINGW_PATH}/bin/"
#             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/nsi_bin/"
#             FILES_MATCHING PATTERN "libstdc++*.dll")
#endif()

# ###############################################################################
# # PACKAGING
# ###############################################################################
# # PACKAGING OPTIONS: GENERAL
set(PACKAGE_NICKNAME "percolator")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Percolator and qvality, two post processors for shotgun proteomics data.")
include("CPack.txt")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxerces-c-dev, libboost-dev, libboost-filesystem-dev (>=1.46), libboost-system-dev")
set(CPACK_RPM_PACKAGE_DEPENDS "libc6, libgcc1")
set(CPACK_RPM_PACKAGE_DEPENDS "${CPACK_RPM_PACKAGE_DEPENDS}, xerces-c-devel, boost-devel, boost-filesystem (>=1.46), boost-system")


include(CPack)
