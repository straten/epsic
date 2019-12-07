//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2016 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/Spinor.h

#ifndef __Spinor_H
#define __Spinor_H

#include "Vector.h"
#include "Jones.h"

#include <iostream>

/***************************************************************************
 *
 *  Spinor class represents an instance of the electric field
 *
 ***************************************************************************/

template<typename T>
class Spinor
{
public:
  Spinor (const std::complex<T>& _x, const std::complex<T>& _y) : x(_x), y(_y){}
  Spinor () { }

  std::complex<T> x;
  std::complex<T> y;

  template <typename U>
  const Spinor& operator *= (U scale) { x *= scale; y *= scale; return *this; }
  const Spinor& operator /= (T norm) { x /= norm; y /= norm; return *this; }
  const Spinor& operator += (const Spinor& e) { x+=e.x; y+=e.y; return *this; }
};

template <typename T>
std::ostream& operator << (std::ostream& os, const Spinor<T>& s)
{
  os << s.x << " " << s.y;
  return os;
}

template <typename T>
const Spinor<T> operator + (Spinor<T> s, const Spinor<T>& t)
{
  return s += t;
}

//! returns the input Spinor transformed by the Jones matrix
template<typename T>
const Spinor<T> operator * (const Jones<T>& j, const Spinor<T>& in)
{
  return Spinor<T> ( j.j00 * in.x + j.j01 * in.y,
		     j.j10 * in.x + j.j11 * in.y );
}

template<typename T, typename U>
const Spinor<T> operator * (U a, Spinor<T> in)
{
  return in *= a;
}

template<typename T, typename U>
const Spinor<T> operator * (Spinor<T> in, U a)
{
  return in *= a;
}


template<typename T, typename U>
void compute_stokes (Vector<4,T>& stokes, const Spinor<U>& e)
{
  double var_x = norm(e.x);
  double var_y = norm(e.y);

  std::complex<double> c_xy = conj(e.x) * e.y;
  
  stokes[0] = var_x + var_y;
  stokes[1] = var_x - var_y;
  stokes[2] = 2.0*c_xy.real();
  stokes[3] = 2.0*c_xy.imag();
}

#endif
