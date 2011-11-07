SET(CMAKE_SYSTEM_NAME Windows)

# specify the cross compiler
find_file(MINGWV NAMES i586-mingw32msvc-gcc i586-mingw32msvc-g++ PATHS
			/usr/bin
			/usr/local/bin)
if(NOT MINGWV)
  set(MINGW_PREFIX "/usr/bin/i686-pc-mingw32-")
  if(NOT DEFINED ${MING_PATH})
    set(MING_PATH "/usr/i686-pc-mingw32/sys-root/mingw")
  endif()
  SET(Boost_COMPILER "-gcc45")
else(NOT MINGWV)
  set(MINGW_PREFIX "/usr/bin/i586-mingw32msvc-")
  if(NOT DEFINED ${MING_PATH})
    set(MING_PATH "/usr/i586-mingw32msvc")
  endif()
  SET(Boost_COMPILER "-mgw")
endif(NOT MINGWV)

SET(MINGWLIB "${MING_PATH}/lib")

# IF (MSVC)
#   SET (BOOST_LIB_PREFIX "lib")
#   SET (BOOST_COMPILER "-vc80") #default for windows compilation
#   SET(CMAKE_FIND_ROOT_PATH "
#   IF (MSVC71)
#    SET (BOOST_COMPILER "-vc71")
#   ENDIF(MSVC71)
#   IF (MSVC80)
#    SET (BOOST_COMPILER "-vc80")
#   ENDIF(MSVC80)
#   IF (MSVC90)
#    SET (BOOST_COMPILER "-vc90")
#   ENDIF(MSVC90)
#   IF (MSVC91)
#    SET (BOOST_COMPILER "-vc91")
#   ENDIF(MSVC91)
#   IF (MSVC100)
#    SET (BOOST_COMPILER "-vc100")
#   ENDIF(MSVC100)
# ENDIF( MSVC )

# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mwindows -mthreads -static-libgcc")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -mwindows -mthreads -static-libgcc")
# SET(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined,--fix-cortex-a8,-lsupc++ -lstdc++ -L${CMAKE_INSTALL_PREFIX}/lib")
# SET(CMAKE_MODULE_LINKER_FLAGS "-Wl,--no-undefined,--fix-cortex-a8,-lsupc++ -lstdc++ -L${CMAKE_INSTALL_PREFIX}/lib")

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


