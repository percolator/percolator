#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT
AC_CONFIG_SRCDIR([main.cpp])
AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE(percolator, 1.04)

# Debug
# add checking message
AC_MSG_CHECKING(whether to build with debug information)

# create configure parameter and init variable $debuger
debuger=no
AC_ARG_ENABLE(debug,
	AC_HELP_STRING(
	    [--enable-debug],
	    [Turn on debuging (def=no)]
	),
  	debuger="$enableval"
)

# resulting value to screen (to complete checking message)
AC_MSG_RESULT($debuger)

# set DEBUG flag and introduce additional compiler flags
if test x"$debuger" = x"yes"; then
	AC_DEFINE([DEBUG],[],[Flag indicating that debug is enabled])
	CXXFLAGS = "-g -Wall -Werror"
	CFLAGS   = "-g -Wall -Werror"
fi
	
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
