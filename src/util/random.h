//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2003 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/random.h

#ifndef __random_H
#define __random_H

#include <complex>
#include <inttypes.h>

// returns the current microsecond
uint64_t usec_seed ();

// seeds the random number generator with the current microsecond
void random_init ();

// uniformly distributed on 0,1
double random_double ();

template <class T, class U>
void random_value (T& value, U scale)
{
  value = (random_double() - 0.5) * 2.0 * scale;
}

template <class T, class U>
void random_value (std::complex<T>& value, U scale)
{
  T real=0, imag=0;

  random_value (real, scale);
  random_value (imag, scale);

  value = std::complex<T> (real, imag);
}

template <class T, class U>
void random_vector (T& array, U scale)
{
  for (unsigned i=0; i<array.size(); i++)
    random_value (array[i], scale);
}

template <class T, class U>
void random_matrix (T& array, U scale)
{
  for (unsigned i=0; i<array.size(); i++)
    random_vector (array[i], scale);
}

#endif
