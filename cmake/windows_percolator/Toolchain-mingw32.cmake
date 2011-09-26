SET(CMAKE_SYSTEM_NAME Windows)

# specify the cross compiler
# SET(CMAKE_C_COMPILER /usr/bin/i686-pc-mingw32-gcc)
# SET(CMAKE_CXX_COMPILER /usr/bin/i686-pc-mingw32-g++)
find_file(MINGWV NAMES i586-mingw32msvc-gcc i586-mingw32msvc-g++ PATHS
			/usr/bin
			/usr/local/bin)
if(MINGWV-NOTFOUND)
  set(MINGW_PREFIX "i686-pc-mingw32-")
else(MINGWV-NOTFOUND)
  set(MINGW_PREFIX "i586-mingw32msvc-")
endif()
set(CMAKE_AR ${MINGW_PREFIX}ar)
set(CMAKE_RANLIB ${MINGW_PREFIX}ranlib)
set(CMAKE_LINKER ${MINGW_PREFIX}ld)
set(CMAKE_C_COMPILER ${MINGW_PREFIX}gcc)
set(CMAKE_CXX_COMPILER ${MINGW_PREFIX}g++)
set(CMAKE_EXECUTABLE_SUFFIX ".exe")
# where is the target environment
SET(CMAKE_FIND_ROOT_PATH /usr/i686-pc-mingw32;/usr/i586-mingw32msvc;/usr/i686-pc-mingw32/sys-root/mingw;/usr/i586-mingw32msvc/sys-root/mingw;
			  ${PROJECT_SOURCE_DIR}/../libs;${PROJECT_SOURCE_DIR}/libs;${PROJECT_SOURCE_DIR}/../../libs;
			  ${PROJECT_SOURCE_DIR}/libs/xsd-3.0.0-i686-windows/libxsd;${PROJECT_SOURCE_DIR}/../libs/xsd-3.0.0-i686-windows/libxsd
			  ${PROJECT_SOURCE_DIR}/../../libs/xsd-3.0.0-i686-windows/libxsd)
# search for programs in the build host directories
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM BOTH)
# for libraries and headers in the target directories
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

# FindQt4.cmake queries qmake to get information,
# which doesn't work when crosscompiling
# SET(QT_HEADERS_DIR ${CMAKE_FIND_ROOT_PATH}/include)
# SET(QT_LIBRARY_DIR ${CMAKE_FIND_ROOT_PATH}/lib)

# #----------------------------------------------------------
# # disable dynamic links
# #----------------------------------------------------------
# SET(CMAKE_SKIP_RPATH TRUE)
# SET(CMAKE_SHARED_LIBRARY_C_FLAGS "")            # -pic
# SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "")       # -shared
# SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")         # +s, flag for exe link to use shared lib
# SET(CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG "")       # -rpath
# SET(CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG_SEP "")   # : or empty
# 
# SET(CMAKE_LINK_LIBRARY_SUFFIX "")
# SET(CMAKE_STATIC_LIBRARY_PREFIX "lib")
# SET(CMAKE_STATIC_LIBRARY_SUFFIX ".a")
# SET(CMAKE_SHARED_LIBRARY_PREFIX "lib")          # lib
# SET(CMAKE_SHARED_LIBRARY_SUFFIX ".a")           # .a
# SET(CMAKE_EXECUTABLE_SUFFIX "")          # .exe
# SET(CMAKE_DL_LIBS "" )
# 
# SET(CMAKE_FIND_LIBRARY_PREFIXES "lib")
# SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# 
# SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
# SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
# SET(CMAKE_C_LINK_SHARED_LIBRARY)
# SET(CMAKE_C_LINK_MODULE_LIBRARY)
