//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2025 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/epsic/complex.h

#ifndef _epsic_complex_type_H
#define _epsic_complex_type_H

#ifdef __CUDACC__
#include <cuda/std/complex>
#else
#include <complex>
#endif

namespace epsic
{
#ifdef __CUDACC__
  template <typename T>
  using complex = cuda::std::complex<T>;
#else
  template <typename T>
  using complex = std::complex<T>;
#endif
}

#ifdef __CUDACC__
#define __prefix__ __host__ __device__
#else
#define __prefix__
#endif

#endif

