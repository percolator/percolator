#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

echo "Combined configure script for qvality and percolator"

AC_PREREQ(2.59)
AC_INIT
AC_CONFIG_SRCDIR([qvality_main.cpp])
AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER([config.h])

# enable percolator installation
AC_MSG_CHECKING(whether to build percolator)

buildpercolator=no
if test -e "${srcdir}/Caller.cpp"
then
  AC_ARG_ENABLE(percolator,
      AC_HELP_STRING(
      [--enable-percolator],
      [Build percolator applications. Default is enabled.]
      ),
      buildpercolator="$enableval", buildpercolator=yes)
fi

AC_MSG_RESULT($buildpercolator)
AM_CONDITIONAL(PERCOLATOR, test [x"$buildpercolator" = x"yes"])

if test x"$buildpercolator" = x"yes"; then
  mypack="percolator"
else
  mypack="qvality"
fi

AM_INIT_AUTOMAKE($mypack, "xiaoyu2")

## Debug
# add checking message
# save $CFLAGS since AC_PROG_CC likes to insert "-g -O2"
# if $CFLAGS is blank
cflags_save="$CFLAGS"
cppflags_save="$CPPFLAGS"
cxxflags_save="$CXXFLAGS"
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
CFLAGS="$cflags_save"
CPPFLAGS="$cppflags_save"
CXXFLAGS="$cxxflags_save"

# create configure parameter and init variable $debuger
AC_MSG_CHECKING(whether to build with debug information)
debuger=no
AC_ARG_ENABLE(debug,
	AC_HELP_STRING(
	    [--enable-debug],
	    [Turn on debuging (def=no)]
	),
  	debuger="$enableval")
# resulting value to screen (to complete checking message)
AC_MSG_RESULT($debuger)

# create configure parameter and init variable $boundcheck
AC_MSG_CHECKING(whether to build with bound checks for the gcc stl containers)
boundscheck=no
AC_ARG_ENABLE(boundscheck,
	AC_HELP_STRING(
	    [--enable-boundscheck],
	    [Turn on boundscheck for GCC STL library (def=no)]
	),
  	boundscheck="$enableval")
# resulting value to screen (to complete checking message)
AC_MSG_RESULT($boundscheck)


# set DEBUG flag and introduce additional compiler flags
if test x"$debuger" = x"yes"; then
    DEBUG_OPTIONS="-ggdb -O0 -Wall -Werror -DDEBUG"
	AC_DEFINE([DEBUG],[],[Flag indicating that debug options are enabled])
else
    DEBUG_OPTIONS="-O3"
fi

if test x"$boundscheck" = x"yes"; then
    DEBUG_OPTIONS="$DEBUG_OPTIONS -D_GLIBCXX_DEBUG"
	AC_DEFINE([_GLIBCXX_DEBUG],[],[Flag indicating that GCC STL container bound check is enabled])
fi

AC_SUBST([DEBUG_OPTIONS])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([sqrt])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
