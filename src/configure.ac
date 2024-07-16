#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)

# AC_INIT version should match version-info in Makefile.am
AC_INIT([EPSIC],[2.0.0],[vanstraten.willem@gmail.com])

AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_SRCDIR([epsic.cpp])

# Initialize automake to place object files in sub-directories
AM_INIT_AUTOMAKE([foreign subdir-objects])

AC_PROG_CXX
AC_PROG_CXXCPP

EPSIC_BEST_PARTIAL_SPECIALIZATION
EPSIC_DEFAULT_PARTIAL_SPECIALIZATION
EPSIC_COMPLEX_TEMPLATES

LT_INIT
AC_SUBST([LIBTOOL_DEPS])

AC_PROG_INSTALL

SWIN_LIB_HEALPIX

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([util/Makefile Makefile])

AC_OUTPUT

echo
echo "EPSIC is now ready to be compiled.  Please run 'make'"
echo

