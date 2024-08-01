//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2004 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/complex_math.h

#ifndef _epsic_complex_math_H
#define _epsic_complex_math_H

#include <complex>
#include "true_math.h"

namespace true_math
{
  template<typename T>
  bool finite (const std::complex<T>& z)
  { return true_math::finite(z.real()) && true_math::finite(z.imag()); }
}

#endif

