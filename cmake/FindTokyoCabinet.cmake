# This module defines
#
#  TOKYOCABINET_FOUND - Tokyo Cabinet was found
#  TOKYOCABINET_INCLUDE_DIR - The Tokyo Cabinet include directory
#  TOKYOCABINET_LIBRARIES - The Tokyo Cabinet libraries


find_path(TOKYOCABINET_INCLUDE_DIR tcbdb.h )
find_library(TOKYOCABINET_LIBRARIES NAMES tokyocabinet )
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TokyoCabinet DEFAULT_MSG TOKYOCABINET_LIBRARIES TOKYOCABINET_INCLUDE_DIR)


