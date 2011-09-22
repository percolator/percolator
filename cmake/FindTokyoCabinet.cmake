# This module defines
#
#  TOKYOCABINET_FOUND - Tokyo Cabinet was found
#  TOKYOCABINET_INCLUDE_DIR - The Tokyo Cabinet include directory
#  TOKYOCABINET_LIBRARIES - The Tokyo Cabinet libraries


find_path(TOKYOCABINET_INCLUDE_DIR tcutil.h tcbdb.h )
find_library(TOKYOCABINET_LIBRARIES NAMES tokyocabinet libtokyocabinet PATH )

IF(TOKYOCABINET_INCLUDE_DIR AND TOKYOCABINET_LIBRARY)
  SET(TOKYOCABINET_FOUND TRUE)
ENDIF(TOKYOCABINET_INCLUDE_DIR AND TOKYOCABINET_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TokyoCabinet DEFAULT_MSG TOKYOCABINET_LIBRARIES TOKYOCABINET_INCLUDE_DIR)


