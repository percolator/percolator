SET(CMAKE_SYSTEM_NAME Windows)

# specify the cross compiler
find_file(MINGWV NAMES i586-mingw32msvc-gcc i586-mingw32msvc-g++ PATHS
			/usr/bin
			/usr/local/bin)
if(NOT MINGWV)
  set(MINGW_PREFIX "/usr/bin/i686-pc-mingw32-")
  if(NOT DEFINED ${MING_PATH})
    set(MINGW_PATH "/usr/i686-pc-mingw32/sys-root/mingw")
  endif()
  if(NOT DEFINED ${Boost_COMPILER})
    SET(Boost_COMPILER "-gcc45")
  endif()
else(NOT MINGWV)
  set(MINGW_PREFIX "/usr/bin/i586-mingw32msvc-")
  if(NOT DEFINED ${MING_PATH})
    set(MINGW_PATH "/usr/i586-mingw32msvc")
  endif()
  if(NOT DEFINED ${Boost_COMPILER})
    SET(Boost_COMPILER "-mgw")
  endif()
endif(NOT MINGWV)

SET (BOOST_LIB_PREFIX "lib")
##probably not necessary, need to set MINGWLIB accordingly and add the linker directory
set(MINGWLIB "${MING_PATH}/lib")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -static-libgcc")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -static-libgcc")

set(CMAKE_C_COMPILER    ${MINGW_PREFIX}gcc)
set(CMAKE_CXX_COMPILER  ${MINGW_PREFIX}g++)
set(CMAKE_RC_COMPILER   ${MINGW_PREFIX}windres)
# set(CMAKE_AR 		${MINGW_PREFIX}ar)
set(CMAKE_RANLIB 	${MINGW_PREFIX}ranlib)
set(CMAKE_LINKER 	${MINGW_PREFIX}ld)
set(CMAKE_STRIP         ${MINGW_PREFIX}strip)
set(CMAKE_EXECUTABLE_SUFFIX ".exe")

# Boost Configuration
SET(Boost_USE_STATIC_LIBS ON) 

# where is the target environment
SET(CMAKE_FIND_ROOT_PATH /usr/i686-pc-mingw32;/usr/i586-mingw32msvc;/usr/i686-pc-mingw32/sys-root/mingw;/usr/i586-mingw32msvc/sys-root/mingw;
			 ${MINGW_PATH};${MINGW_PATH}/sys-root/mingw;
			 $ENV{HOME}/i586-mingw32msvc;$ENV{HOME}/i686-pc-mingw32;${CMAKE_PREFIX_PATH};
			 ${PROJECT_SOURCE_DIR}/../libs;${PROJECT_SOURCE_DIR}/libs;${PROJECT_SOURCE_DIR}/../../libs)

# search for programs in the build host directories
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
# for libraries and headers in the target directories
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)


