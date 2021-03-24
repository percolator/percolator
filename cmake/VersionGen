# Adapted from https://github.com/yasm/yasm/blob/master/cmake/modules/VersionGen.cmake (2016-12-16)
# Edited for use in Percolator by MT, 2016
#

macro (TAG_TO_VERSION _nightly _version _tag)
  string (REPLACE "rel-" "" _vn_cmake_format "${_tag}")
  string (REPLACE "-" "." _vn_cmake_format "${_vn_cmake_format}")
  set(_version_tmp ${_vn_cmake_format})
  set(_nightly_tmp "0")
  if (_vn_cmake_format MATCHES "^[0-9]*\\.[0-9]*\\.[0-9]*\\.g.*") # patch number is not present in tag description
    string( REGEX REPLACE "([0-9]*)\\.([0-9]*)\\.[0-9]*\\.g.*" "\\1.\\2.0" _version_tmp ${_vn_cmake_format} )
    string( REGEX REPLACE "[0-9]*\\.[0-9]*\\.([0-9]*)\\.g(.*)" "nightly-\\1-\\2" _nightly_tmp ${_vn_cmake_format} )
  elseif (_vn_cmake_format MATCHES "^[0-9]*\\.[0-9]*\\.[0-9]*\\.[0-9]*\\.g.*") # patch number is present in tag description
    string( REGEX REPLACE "([0-9]*)\\.([0-9]*)\\.([0-9]*)\\.[0-9]*\\.g.*" "\\1.\\2.\\3" _version_tmp ${_vn_cmake_format} )
    string( REGEX REPLACE "[0-9]*\\.[0-9]*\\.[0-9]*\\.([0-9]*)\\.g(.*)" "nightly-\\1-\\2" _nightly_tmp ${_vn_cmake_format} )
  elseif (_vn_cmake_format MATCHES "^[0-9]*\\.[0-9]*\\.[0-9]*$") # patch release
    string( REGEX REPLACE "[0-9]*\\.[0-9]*\\.([0-9]*)" "\\1" _nightly_tmp ${_vn_cmake_format} )
  endif (_vn_cmake_format MATCHES "^[0-9]*\\.[0-9]*\\.[0-9]*\\.g.*")
  set(${_version} ${_version_tmp})
  set(${_nightly} ${_nightly_tmp})
endmacro (TAG_TO_VERSION)

macro (NIGHTLY_VERSION_GEN _nightly_version _major_version _minor_version _patch_version)
  set (_vn "${_patch_version}")
  
  set(_release_tag "rel-${_major_version}-${_minor_version}")
  if (NOT "${_patch_version}" STREQUAL "0")
    set(_release_tag "${_release_tag}-${_patch_version}")
  endif (NOT "${_patch_version}" STREQUAL "0")

  if (EXISTS "${PERCOLATOR_SOURCE_DIR}/.git")
    execute_process (COMMAND git describe --tags --match "rel-[0-9]*" HEAD
                     RESULT_VARIABLE _git_result
                     OUTPUT_VARIABLE _git_vn
                     ERROR_QUIET
                     WORKING_DIRECTORY ${PERCOLATOR_SOURCE_DIR}
                     OUTPUT_STRIP_TRAILING_WHITESPACE)
    #message(STATUS "${_git_result} ${_git_vn} ${_release_tag}")
    
    if (_git_result EQUAL 0)
        
      # extract version number string from tag descriptions
      TAG_TO_VERSION(_tmp _current_version "${_release_tag}")
      TAG_TO_VERSION(_nightly_tag _last_tagged_version "${_git_vn}")
      #message(STATUS "${_current_version}, ${_tmp}, ${_last_tagged_version}, ${_nightly_tag}")
      
      # if we are trying to build a new release, do not use the nightly and/or dirty suffix
      if (NOT "${_current_version}" VERSION_GREATER "${_last_tagged_version}")
        set(_vn ${_nightly_tag})
        
        # Append -dirty if there are local changes
        execute_process (COMMAND git update-index -q --refresh
                         WORKING_DIRECTORY ${PERCOLATOR_SOURCE_DIR})
        execute_process (COMMAND git diff-index --name-only HEAD --
                         OUTPUT_VARIABLE _git_vn_dirty
                         WORKING_DIRECTORY ${PERCOLATOR_SOURCE_DIR})
        if (_git_vn_dirty)
            set (_vn "${_vn}-dirty")
        endif (_git_vn_dirty)
          
      endif(NOT "${_current_version}" VERSION_GREATER "${_last_tagged_version}")

    endif (_git_result EQUAL 0)
  endif (EXISTS "${PERCOLATOR_SOURCE_DIR}/.git")

  # Set output version variable
  set (${_nightly_version} ${_vn})
endmacro (NIGHTLY_VERSION_GEN)
