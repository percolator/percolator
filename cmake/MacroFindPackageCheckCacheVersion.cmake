# MACRO_FIND_PACKAGE_WITH_CACHE_AND_VERSION_CHECK(pkg name header [lib ...])
#
# The same with MACRO_FIND_PACKAGE_WITH_CACHE. Version is also checked. If the
# version does not match, the package is considered not found.
#
# User is responsibile to provide a function GET_PACKAGE_VERSION. The function
# can use variables such as ${name}_INCLUDE_DIR and returns the version string.

MACRO(MACRO_FIND_PACKAGE_WITH_CACHE_AND_VERSION_CHECK pkg name header)

  SET(_${name}_IN_CACHE FALSE)

  IF(${name}_INCLUDE_DIR)
    SET(_${name}_IN_CACHE TRUE)
    FOREACH(_lib ${ARGN})
      IF(NOT ${name}_${_lib}_LIBRARY)
        SET(_${name}_IN_CACHE FALSE)
      ENDIF(NOT ${name}_${_lib}_LIBRARY)
    ENDFOREACH(_lib)

    IF(_${name}_IN_CACHE)
      IF(NOT ${name}_VERSION)
        GET_PACKAGE_VERSION()
      ENDIF(NOT ${name}_VERSION)

      INCLUDE(MacroFindPackageCheckCacheVersion)
      MACRO_FIND_PACKAGE_CHECK_CACHE_VERSION(_${name}_IN_CACHE ${name})
    ENDIF(_${name}_IN_CACHE)

    IF(NOT _${name}_IN_CACHE)
      SET(${name}_INCLUDE_DIR) # remove it
    ENDIF(NOT _${name}_IN_CACHE)
  ENDIF(${name}_INCLUDE_DIR)

  IF(_${name}_IN_CACHE)
    # in cache already
    MESSAGE(STATUS "${pkg} in cache")

    SET(${name}_FOUND TRUE)

    SET(${name}_INCLUDE_DIRS ${${name}_INCLUDE_DIR})
    SET(${name}_LIBRARIES ${${name}_LIBRARY})
    FOREACH(_lib ${ARGN})
      IF(${name}_${_lib}_LIBRARY)
        LIST(APPEND ${name}_LIBRARIES ${${name}_${_lib}_LIBRARY})
      ENDIF(${name}_${_lib}_LIBRARY)
    ENDFOREACH(_lib)

  ELSE(_${name}_IN_CACHE)
    # should search it

    INCLUDE(MacroFindPackage)
    MACRO_FIND_PACKAGE(${pkg} ${name} ${header} ${ARGN})

    IF(${name}_FOUND)
      GET_PACKAGE_VERSION()
      IF(NOT ${name}_VERSION)
        SET(${name}_FOUND FALSE)
      ENDIF(NOT ${name}_VERSION)
    ENDIF(${name}_FOUND)

    # checks version if user specified one
    SET(_details "Find ${pkg}: failed.")
    IF(${name}_FOUND AND ${name}_FIND_VERSION)
      INCLUDE(MacroVersionCmp)
      MACRO_VERSION_CMP("${${name}_VERSION}" "${${name}_FIND_VERSION}" _cmp_result)
      IF(_cmp_result LESS 0)
        SET(_details "${_details} ${${name}_FIND_VERSION} required but ${${name}_VERSION} found")
        SET(${name}_FOUND FALSE)
      ELSEIF(${name}_FIND_VERSION_EXACT AND _cmp_result GREATER 0)
        SET(_details "${_details} exact ${${name}_FIND_VERSION} required but ${${name}_VERSION} found")
        SET(${name}_FOUND FALSE)
      ENDIF(_cmp_result LESS 0)
    ENDIF(${name}_FOUND AND ${name}_FIND_VERSION)

    IF(NOT ${name}_FOUND)
      IF(${name}_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "${_details}")
      ELSEIF(NOT ${name}_FIND_QUIETLY)
        MESSAGE(STATUS "${_details}")
      ENDIF(${name}_FIND_REQUIRED)
    ENDIF(NOT ${name}_FOUND)

  ENDIF(_${name}_IN_CACHE)

ENDMACRO(MACRO_FIND_PACKAGE_WITH_CACHE_AND_VERSION_CHECK)