# Redistribution and use is allowed according to the terms of the BSD license.
#
# Copyright (c) 2011 Peter Johnson
#
# Adapted from https://github.com/yasm/yasm/blob/master/cmake/modules/VersionGen.cmake (2016-12-16)
# Edited for use in Percolator by MT, 2016
#

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
        message(STATUS "${_git_result} ${_git_vn} ${_release_tag}")
        
        if (_git_result EQUAL 0)
        
            if (NOT "${_release_tag}" STREQUAL "${_git_vn}")
                # Special handling until we get a more recent tag on the master
                # branch
                if (_git_vn MATCHES "^rel-[0-9-]*-g")
                    execute_process (COMMAND git merge-base ${_release_tag} HEAD
                                     OUTPUT_VARIABLE _merge_base
                                     ERROR_QUIET
                                     WORKING_DIRECTORY ${PERCOLATOR_SOURCE_DIR}
                                     OUTPUT_STRIP_TRAILING_WHITESPACE)
                    message (STATUS "Merge base: ${_merge_base}")

                    execute_process (COMMAND git rev-list ${_merge_base}..HEAD
                                     OUTPUT_VARIABLE _rev_list
                                     WORKING_DIRECTORY ${PERCOLATOR_SOURCE_DIR}
                                     ERROR_QUIET)
                    string (REGEX MATCHALL "[^\n]*\n" _rev_list_lines "${_rev_list}")
                    #message (STATUS "Rev list: ${_rev_list_lines}")
                    list (LENGTH _rev_list_lines _vn1)

                    execute_process (COMMAND git rev-list --max-count=1 --abbrev-commit HEAD
                                     OUTPUT_VARIABLE _vn2
                                     ERROR_QUIET
                                     WORKING_DIRECTORY ${PERCOLATOR_SOURCE_DIR}
                                     OUTPUT_STRIP_TRAILING_WHITESPACE)

                    set (_vn "nightly-commit${_vn1}-${_vn2}")
                endif (_git_vn MATCHES "^rel-[0-9-]*-g")
                
            endif(NOT "${_release_tag}" STREQUAL "${_git_vn}")

            # Append -dirty if there are local changes
            execute_process (COMMAND git update-index -q --refresh
                             WORKING_DIRECTORY ${PERCOLATOR_SOURCE_DIR})
            execute_process (COMMAND git diff-index --name-only HEAD --
                             OUTPUT_VARIABLE _git_vn_dirty
                             WORKING_DIRECTORY ${PERCOLATOR_SOURCE_DIR})
            if (_git_vn_dirty)
                set (_vn "${_vn}-dirty")
            endif (_git_vn_dirty)

        endif (_git_result EQUAL 0)
    endif (EXISTS "${PERCOLATOR_SOURCE_DIR}/.git")

    # Set output version variable
    set (${_nightly_version} ${_vn})
endmacro (NIGHTLY_VERSION_GEN)
