#
# SWIN_LIB_HEALPIX([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
#
# This m4 macro checks availability of the HEALPIX Library
#
# HEALPIX_CFLAGS - autoconfig variable with flags required for compiling
# HEALPIX_LIBS   - autoconfig variable with flags required for linking
# HAVE_HEALPIX   - automake conditional
# HAVE_HEALPIX   - pre-processor macro in config.h
#
# Notice that the environment variable HEALPIX is required
#
# ----------------------------------------------------------
AC_DEFUN([SWIN_LIB_HEALPIX],
[
  AC_PROVIDE([SWIN_LIB_HEALPIX])

  AC_MSG_CHECKING([for HEALPix libary installation])

  HEALPIX_CFLAGS="-I${HEALPIX}/include/healpix_cxx"
  HEALPIX_LIBS="-L${HEALPIX}/lib -lhealpix_cxx -lsharp"

  ac_save_CXXFLAGS="$CXXFLAGS"
  ac_save_LIBS="$LIBS"

  AC_LANG_PUSH(C++)

  LIBS="$ac_save_LIBS $HEALPIX_LIBS"
  CXXFLAGS="$ac_save_CXXFLAGS $HEALPIX_CFLAGS"

  AC_TRY_LINK([#include "healpix_base.h"
               #include "healpix_map.h"],
              [Healpix_Map<double> map = Healpix_Map<double>(); 
               map.SetNside(128, RING); ],
              have_healpix=yes, have_healpix=no)

  AC_MSG_RESULT($have_healpix)

  LIBS="$ac_save_LIBS"
  CXXFLAGS="$ac_save_CXXFLAGS"

  AC_LANG_POP(C++)

  if test x"$have_healpix" = xyes; then
    AC_DEFINE([HAVE_HEALPIX], [1], [Define to 1 if you have the HEALPix library])
    [$1]
  else
    AC_MSG_WARN([HEALPix code will not be compiled])
    if test x"$HEALPIX" = x; then
      AC_MSG_WARN([Please set the HEALPIX environment variable])
    fi
    HEALPIX_CFLAGS=""
    HEALPIX_LIBS=""
    [$2]
  fi

  AC_SUBST(HEALPIX_CFLAGS)
  AC_SUBST(HEALPIX_LIBS)
  AM_CONDITIONAL(HAVE_HEALPIX, [test x"$have_healpix" = xyes])

])

