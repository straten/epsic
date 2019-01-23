#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([EPSIC],[2019-01-23],[vanstraten.willem@gmail.com])

AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_SRCDIR([Util/Makefile.am])

AC_ARG_VAR([PSRHOME],    [Standard pulsar home directory])
AC_ARG_VAR([LOGIN_ARCH], [Architecture-dependent sub-directory of PSRHOME])
AC_ARG_VAR([PACKAGES],   [Root directory of third-party packages])
AC_ARG_VAR([PGPLOT_DIR], [Pgplot graphics library directory])

AC_PREFIX_DEFAULT([${PSRHOME:-"/usr/local"}${PSRHOME:+"/$LOGIN_ARCH"}])

################################################################################
#
# Warn about pulsar group environment variables
#

if test -z "$PSRHOME"; then
  AC_MSG_WARN([The PSRHOME environment variable is not set])
else
  if test -z "$LOGIN_ARCH"; then
    AC_MSG_WARN([The LOGIN_ARCH environment variable is not set])
  fi
fi


################################################################################
#
# Set the @prefix@ variable in autoconf-generated files
#
if test "x$prefix" = xNONE; then
  prefix=$ac_default_prefix
fi

PSRCHIVE_INSTALL=$prefix
PSRCHIVE_HTTP="http://psrchive.sourceforge.net"

AC_SUBST(PSRCHIVE_INSTALL)
AC_SUBST(PSRCHIVE_HTTP)

# Check if the user has set compiler options
SWIN_OPTIONS_SET

# Enable convenience compiler selection
SWIN_COMPILER

# Initialize automake to place object files in sub-directories
AM_INIT_AUTOMAKE([foreign subdir-objects])

# Create static libraries by default
AC_DISABLE_SHARED
AC_SUBST(enable_shared)

# Ensure that SIGPROC find does not over-ride the standard
SWIN_PROG_FIND

# Checks for programs.

AC_PROG_CXX
SWIN_BEST_F77

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXXCPP
AM_PROG_CC_C_O

AC_PROG_INSTALL
AC_PROG_LIBTOOL

AC_F77_WRAPPERS
AC_F77_LIBRARY_LDFLAGS

# Determine the machine endian
AC_C_BIGENDIAN

# Enable large-file support
AC_SYS_LARGEFILE

# Check which version of the pulsar catalog is installed
SWIN_PSRCAT

#
# Disable the debugging information option, -g
#
SWIN_DEBUG
SWIN_LOCAL

ATNF_DISABLE_OPTIMIZE

#
# Check C++ template compilation
#
SWIN_BEST_PARTIAL_SPECIALIZATION
SWIN_DEFAULT_PARTIAL_SPECIALIZATION
SWIN_COMPLEX_TEMPLATE
SWIN_TEST_LIMITS

#
# Set up to use any packages already installed in the target directory
#
AC_ARG_ENABLE([existing],
              AC_HELP_STRING([--disable-existing],
              [Don't look for dependencies in existing installation]))

if test -d $PSRCHIVE_INSTALL -a x"$enable_existing" != x"no"; then
  CPPFLAGS="-I$PSRCHIVE_INSTALL/include $CPPFLAGS"
  LDFLAGS="-L$PSRCHIVE_INSTALL/lib $LDFLAGS"
fi

#
# enable POSIX C99 format macros
#
AC_DEFINE([__STDC_FORMAT_MACROS], [1], [enable POSIX C99 format macros])

#
# Checks for header files.
#
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([malloc.h sys/statvfs.h sys/vfs.h sys/mount.h])

#
# Checks for essential libraries
#
AC_CHECK_LIB([m], [sin])
SWIN_LIB_FFTW
SWIN_LIB_MKL
SWIN_LIB_IPP

#
# Set up for automatic package installation
#
SWIN_DOWNLOAD
AC_CONFIG_FILES([packages/Makefile
	packages/cfitsio.csh
	packages/fftw.csh
	packages/pgplot.csh
	packages/psrcat.csh
	packages/tempo.csh
	packages/tempo2.csh])

AC_CONFIG_COMMANDS(executable,[chmod a+x packages/*.csh])

#
# Check that at least one FFT library is available
#

if test x"$have_fftw3" = xno -a x"$have_mkl" = xno -a x"$have_ipp" = xno; then

  # produce the package installation scripts before aborting.
  AC_OUTPUT

  AC_MSG_ERROR([No suitable FFT library was found.

  Please see http://psrchive.sourceforge.net/third/fftw

  or simply run ./packages/fftw.csh to install FFTW automatically
])
fi

AC_SUBST(FFT_LIBS)
FFT_LIBS="$FFTW_LIBS $MKL_LIBS $IPP_LIBS"
AC_SUBST(FFT_CFLAGS)
FFT_CFLAGS="$FFTW_CFLAGS $MKL_CFLAGS $IPP_CFLAGS"

#
# Checks for optional libraries
#
SWIN_LIB_CFITSIO
SWIN_LIB_CULA
SWIN_LIB_PUMA
SWIN_LIB_TEMPO2
SWIN_LIB_MPI
SWIN_LIB_READLINE
SWIN_LIB_GSL
SWIN_LIB_HEALPIX
SWIN_LIB_EIGEN

# Check for psrxml io library
MJK_LIB_PSRXML

#
# Generate python module if --enable-shared is used but not
# if --disable-python is also used.
#
AC_ARG_ENABLE([python],
              AC_HELP_STRING([--disable-python],
              [Don't build python extensions]))
if test x"$enable_shared" == xyes -a x"$enable_python" != xno; then
  AM_PATH_PYTHON(,, [:])
  AC_PROG_SWIG
  SWIG_ENABLE_CXX
  SWIG_PYTHON
else
  PYTHON=":"
fi

AM_CONDITIONAL([HAVE_PYTHON], [test "$PYTHON" != :])

ACX_PTHREAD
LIBS="$PTHREAD_LIBS $LIBS"
CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS"

#
# Compile plotting and GUI code if --enable-lean is not used
#
if test x"$enable_lean" == xyes; then

AM_CONDITIONAL(HAVE_PGPLOT, [test 1 = 0])
AM_CONDITIONAL(HAVE_QT, [test 1 = 0])
AM_CONDITIONAL(HAVE_QTDRIV, [test 1 = 0])

else

#
# Checks for graphics libraries
#
SWIN_LIB_PGPLOT
BNV_HAVE_QT
SWIN_LIB_QTDRIV

fi

#
# For developers: preserve file modification times
#
INSTALL_DATA="${INSTALL} -m 644 -p"
install_sh="CPPROG='cp -p' ${install_sh}"

#
# For developers: provide the ability to compile a lean version of the code
#
AC_ARG_WITH([lean], AC_HELP_STRING([--enable-lean],
				   [compile a lean version of the library]))

AM_CONDITIONAL(MAKE_LEAN, [test x"$enable_lean" == xyes])

#
# Checks for library functions.
#
SWIN_FUNC_GETOPT_LONG


################################################################################
#
#
#

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
		config/Makefile
		Util/Makefile
		Util/third/Makefile
		Util/third/parsifal/Makefile
		Util/third/star/Makefile
		Util/third/gnu/Makefile
		Util/units/Makefile
		Util/fft/Makefile
		Util/fft/install_bench
		Util/genutil/Makefile
		Util/qtutil/Makefile
		Util/fitsutil/Makefile
		Util/tempo/Makefile
		Util/tempo2/Makefile
		Util/resources/Makefile
		Base/Makefile
		Base/Classes/Makefile
		Base/Extensions/Makefile
		Base/Checks/Makefile
		Base/Resources/Makefile
		Base/Formats/Makefile
		Base/Formats/EPN/Makefile
		Base/Formats/EPOS/Makefile
		Base/Formats/Timer/Makefile
		Base/Formats/PSRFITS/Makefile
		Base/Formats/PuMa/Makefile
                Base/Formats/ASCII/Makefile
                Base/Formats/BPP/Makefile
                Base/Formats/WAPP/Makefile
                Base/Formats/PRESTO/Makefile
                Base/Formats/ASP/Makefile
		Base/Formats/UVM/Makefile
		Base/Formats/UVM/uvmlib/Makefile
		Base/Formats/pdv_KtZ/Makefile
		Base/Applications/Makefile
		More/Makefile
		More/General/Makefile
		More/Plotting/Makefile
		More/Plotting/pgutil/Makefile
		More/MEAL/Makefile
		More/Timing/Makefile
		More/Polarimetry/Makefile
		More/Polarimetry/iri2001/Makefile
		More/RFIMitigation/Makefile
		More/Applications/Makefile
		More/Applications/njkk08/Makefile
		More/Simulation/Makefile
		More/rhythm/Makefile
		More/psrgui/Makefile
		More/python/Makefile
		Management/Makefile
		Management/psrchive
		Management/psrchive_install.h
		Management/release.csh
		Makefile])

AC_OUTPUT

if test x"$enable_lean" == xyes; then

   echo
   echo "********************************************"
   echo "A lean version of PSRCHIVE will be compiled."
   echo "This means that the polarimetric calibration"
   echo "software, related error analysis routines,"
   echo "and graphical interfaces will be unavailable."
   echo "********************************************"

else

if test x"$have_tempo2" != xyes; then
   echo
   echo "The TEMPO2 Predictor library was not found: (desirable)"
   echo "  PSRCHIVE applications will not be able to use TEMPO2 predictors."
   echo "  Run ./packages/tempo2.csh to automatically install."
fi

if test x"$have_cfitsio" != xyes; then
   echo
   echo "The CFITSIO library was not found: (desirable)"
   echo "  PSRCHIVE applications will not be able to load PSRFITS files."
   echo "  Run ./packages/cfitsio.csh to automatically install."
fi

if test x"$have_pgplot" != xyes; then
   echo
   echo "The PGPLOT library was not found: (desirable)"
   echo "  PSRCHIVE applications will not be able to plot."
   echo "  Run ./packages/pgplot.csh to automatically install."
fi

if test x"$have_qt" != xyes; then
   echo
   echo "The Qt library was not found: (optional)"
   echo "  PSRCHIVE graphical user interfaces will not be compiled."
fi

if test x"$have_qtdriv" != xyes; then
   echo
   echo "The qtdriv device handler is not installed: (optional)"
   echo "  rhythm will not be compiled."
fi

fi

if test -n "$swig_version"; then
   swig_major=`echo $swig_version | awk -F. '{print $1}'`
   if test x"$swig_major" == x3; then
      echo
      echo "SWIG version 3.x ($swig_version) found:"
      echo "  WARNING, this may be incompatible with PSRCHIVE."
      echo "  If compile errors occur, please try with SWIG 2.x"
      echo "  or use the --disable-python configure option."
      echo
   fi
fi

echo
echo "PSRCHIVE is now ready to be compiled.  Please run 'make'"
echo

