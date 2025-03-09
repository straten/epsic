# EPSIC_LIB_CUDA([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
# ----------------------------------------------------------
AC_DEFUN([EPSIC_LIB_CUDA],
[
  AC_PROVIDE([EPSIC_LIB_CUDA])

  CUDA_CFLAGS=""
  CUDA_LIBS=""

  have_cuda="not found"

  AC_PATH_PROG(cuda_nvcc, nvcc, no)

  AC_MSG_CHECKING([for CUDA installation])

  if test -x $cuda_nvcc; then

    CUDA_BIN=`dirname $cuda_nvcc`
    CUDA_ROOT=`dirname $CUDA_BIN`
    CUDA_CFLAGS=-I$CUDA_ROOT/include
    CUDA_LIBS="-L$CUDA_ROOT/lib64 -lcudart"

    CUDA_NVCC="$cuda_nvcc \$(CUDA_NVCC_FLAGS) -Xcompiler \"\$(DEFAULT_INCLUDES) \$(INCLUDES) \$(AM_CPPFLAGS) \$(CPPFLAGS)\""

    CFLAGS=$CUDA_CFLAGS
    LIBS=$CUDA_LIBS
    AC_TRY_LINK([#include <cuda_runtime.h>], [cudaMalloc (0, 0);],[have_cuda=yes],[have_cuda=no])

  fi

  AC_MSG_RESULT([$have_cuda])

  if test "$have_cuda" = "yes"; then

    AC_DEFINE([HAVE_CUDA],[1],[Define if the CUDA library is present])

  else

    if test "$have_cuda" = "not found"; then
      echo
      AC_MSG_NOTICE([Ensure that CUDA nvcc is in PATH.])
      AC_MSG_NOTICE([Alternatively, use the --with-cuda-dir option.])
      echo
    fi

    [$2]

  fi

  AC_SUBST(CUDA_NVCC)

  CUDA_LIBS="$cuda_LIBS"
  CUDA_CFLAGS="$cuda_CFLAGS"

  AC_SUBST(CUDA_LIBS)
  AC_SUBST(CUDA_CFLAGS)
  AM_CONDITIONAL(HAVE_CUDA,[test "$have_cuda" = "yes"])
])
