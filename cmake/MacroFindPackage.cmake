# - MACRO_FIND_PACKAGE(pkg name header [lib1 ...])
#
# See MACRO_FIND_PACKAGE_WITH_CACHE
# and MACRO_FIND_PACKAGE_WITH_CACHE_AND_VERSION_CHECK
#
# Parameters:
#
# - pkg : package name used for pkg-config
# - name : prefix for result variables
# - header : A header name in the package
# - lib1 ... : list of libraries in the package
#
# The following are set after the configuration is done:
#
# - ${name}_FOUND : Set to True if the package was found
# - ${name}_INCLUDE_DIRS : Include directories
# - ${name}_LIBRARIES : Libraries required to link
#
# User can change following variable in cache to indicate where to find the
# package:
#
# - ${name}_ROOT_DIR : Top directory of the library, typically the header files
# are located in include and libraries are in lib in this
# top directory.
# - ${name}_INCLUDE_DIR,
# ${name}_LIBRARY : Specify include dir and library directly.

MACRO(MACRO_FIND_PACKAGE pkg name header)

  SET(${name}_FOUND FALSE)

  SET(_env_root_dir $ENV{${name}_ROOT_DIR})
  IF(NOT ${name}_ROOT_DIR AND _env_root_dir)
    SET(${name}_ROOT_DIR ${_env_root_dir})
  ENDIF(NOT ${name}_ROOT_DIR AND _env_root_dir)

  # first pkg-config if user does not specify ${name}_INCLUDE_DIR
  IF(NOT ${name}_INCLUDE_DIR AND NOT ${name}_NO_PKG_CONFIG)
    FIND_PACKAGE(PkgConfig)

    IF(PKG_CONFIG_FOUND)
      SET(_backup_pkg_config_path "$ENV{PKG_CONFIG_PATH}")
      IF(${name}_ROOT_DIR)
        SET(ENV{PKG_CONFIG_PATH} "${${name}_ROOT_DIR}/lib/pkgconfig")
      ENDIF(${name}_ROOT_DIR)
      PKG_CHECK_MODULES(${name} ${pkg})
      SET(ENV{PKG_CONFIG_PATH} "${_backup_pkg_config_path}")

      IF(${name}_FOUND)
        SET(${name}_INCLUDE_DIR ${${name}_INCLUDE_DIRS}
          CACHE STRING "${name} header")
        SET(${name}_LIBRARY ${${name}_LIBRARIES}
          CACHE STRING "${name} library")
        MARK_AS_ADVANCED(${name}_LIBRARY)
      ENDIF(${name}_FOUND)

    ENDIF(PKG_CONFIG_FOUND)

  ENDIF(NOT ${name}_INCLUDE_DIR AND NOT ${name}_NO_PKG_CONFIG)

  # then try the normal way
  IF(NOT ${name}_FOUND)

    IF(${name}_ROOT_DIR)
      MESSAGE(STATUS "Try to find ${pkg} in ${${name}_ROOT_DIR}")
      FIND_PATH(${name}_INCLUDE_DIR ${header}
        PATHS "${${name}_ROOT_DIR}/include" "${${name}_ROOT_DIR}/interface"
        NO_DEFAULT_PATH)
      FOREACH(_lib ${ARGN})
        FIND_LIBRARY(${name}_${_lib}_LIBRARY ${_lib}
          PATHS "${${name}_ROOT_DIR}/lib" "${${name}_ROOT_DIR}"
          NO_DEFAULT_PATH)
      ENDFOREACH(_lib)
    ELSE(${name}_ROOT_DIR)
      MESSAGE(STATUS "Try to find ${pkg} in default search locations")
      FIND_PATH(${name}_INCLUDE_DIR ${header})
      FOREACH(_lib ${ARGN})
        FIND_LIBRARY(${name}_${_lib}_LIBRARY ${_lib})
      ENDFOREACH(_lib)
    ENDIF(${name}_ROOT_DIR)

    SET(_lib_list)
    FOREACH(_lib ${ARGN})
      LIST(APPEND _lib_list ${name}_${_lib}_LIBRARY)
      IF(${name}_${_lib}_LIBRARY)
        SET(_lib_found ON)
      ELSE(${name}_${_lib}_LIBRARY)
        SET(_lib_found OFF)
      ENDIF(${name}_${_lib}_LIBRARY)
      SET(${name}_${_lib}_FOUND ${_lib_found}
        CACHE BOOL "Component ${_lib} in ${name} has been found?" FORCE)
      MARK_AS_ADVANCED(${name}_${_lib}_FOUND)
    ENDFOREACH(_lib)
    INCLUDE(FindPackageHandleStandardArgs)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(
      ${name}
      DEFAULT_MSG
      ${name}_INCLUDE_DIR
      ${_lib_list}
      )
    STRING(TOUPPER ${name} _name_upper)
    SET(${name}_FOUND ${${_name_upper}_FOUND})

    IF(${name}_FOUND)
      SET(${name}_LIBRARIES)
      FOREACH(_lib ${_lib_list})
        LIST(APPEND ${name}_LIBRARIES ${${_lib}})
      ENDFOREACH(_lib)
      SET(${name}_INCLUDE_DIRS ${${name}_INCLUDE_DIR})
    ELSE(${name}_FOUND)
      MESSAGE(STATUS "Failed to find ${pkg}")
      MESSAGE(STATUS " Try to set ${name}_ROOT_DIR though environemnt variable or cmake -D")
    ENDIF(${name}_FOUND)

    IF(_lib_list)
      MARK_AS_ADVANCED(${_lib_list})
    ENDIF(_lib_list)

  ENDIF(NOT ${name}_FOUND)

  MARK_AS_ADVANCED(${name}_INCLUDE_DIR)

ENDMACRO(MACRO_FIND_PACKAGE)