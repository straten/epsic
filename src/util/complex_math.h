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

template<typename T>
bool isfinite (const std::complex<T>& z)
{ return isfinite(z.real()) && isfinite(z.imag()); }

#endif

