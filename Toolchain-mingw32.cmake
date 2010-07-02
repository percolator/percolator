SET(CMAKE_SYSTEM_NAME Windows)

# specify the cross compiler
SET(CMAKE_C_COMPILER /usr/bin/i686-pc-mingw32-gcc)
SET(CMAKE_CXX_COMPILER /usr/bin/i686-pc-mingw32-g++)

# where is the target environment
SET(CMAKE_FIND_ROOT_PATH /usr/i686-pc-mingw32/sys-root/mingw)
#SET(CMAKE_FIND_ROOT_PATH /usr/i686-pc-mingw32/sys-root/mingw;/home/esjolund/xsd-3.3.0-i686-windows)
#SET(CMAKE_FIND_ROOT_PATH /home/esjolund/xsd-3.3.0-i686-windows)

# search for programs in the build host directories
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
# for libraries and headers in the target directories
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

# FindQt4.cmake queries qmake to get information,
# which doesn't work when crosscompiling
SET(QT_HEADERS_DIR ${CMAKE_FIND_ROOT_PATH}/include)
SET(QT_LIBRARY_DIR ${CMAKE_FIND_ROOT_PATH}/lib)
