#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([EPSIC],[2019-01-23],[vanstraten.willem@gmail.com])

AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_SRCDIR([epsic.C])

# Initialize automake to place object files in sub-directories
AM_INIT_AUTOMAKE([foreign subdir-objects])

# Create static libraries by default
AC_DISABLE_SHARED
AC_SUBST(enable_shared)

AC_PROG_CXX
AC_PROG_CXXCPP

EPSIC_BEST_PARTIAL_SPECIALIZATION
EPSIC_DEFAULT_PARTIAL_SPECIALIZATION

AC_PROG_INSTALL
AC_PROG_LIBTOOL

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([util/Makefile Makefile])

AC_OUTPUT

echo
echo "EPSIC is now ready to be compiled.  Please run 'make'"
echo

