//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2003 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// espic/src/Vector.h

#ifndef __Vector_H
#define __Vector_H

#include "Traits.h"
#include "Estimate.h"

#include <iostream>
#include <complex>

//! Vector
template <unsigned N, typename T> 
class Vector {

  template<typename U> 
  static void zero (Vector<N, U>& v)
  { for (unsigned i=0; i<N; i++) v.x[i] = T(0.0); }

  //! overloaded template fixes Vector<M,Vector<M,T>> default constructor 
  template<unsigned M, typename U> 
  static void zero (Vector<N, Vector<M, U> >& v)
  { /* do nothing */ }

public:
  T x[N];

  //! Default constructor
  Vector () { zero (*this); }

  Vector (T x0)
  { x[0] = x0; }

  Vector (T x0, T x1)
  { x[0] = x0; x[1] = x1; }

  Vector (T x0, T x1, T x2)
  { x[0] = x0; x[1] = x1; x[2] = x2; }

  Vector (T x0, T x1, T x2, T x3)
  { x[0] = x0; x[1] = x1; x[2] = x2; x[3] = x3; }

  //! Construct from another Vector<U> instance
  template<typename U> Vector (const Vector<N, U>& s)
    { operator=(s); }

  //! Set this instance equal to another Vector<U> instance
  template<typename U> Vector& operator = (const Vector<N, U>& s)
    { for (unsigned i=0; i<N; i++) x[i] = T(s.x[i]); return *this; }

  //! Set this instance equal to a scalar
  Vector& operator = (const T& scalar)
    { x[0] = scalar; for (unsigned i=1; i<N; i++) x[i] = 0.0; return *this; }

  //! Vector addition
  Vector& operator += (const Vector& s)
    { for (unsigned i=0; i<N; i++) x[i] += s.x[i]; return *this; }

  //! Vector subtraction
  Vector& operator -= (const Vector& s)
    { for (unsigned i=0; i<N; i++) x[i] -= s.x[i]; return *this; }

  /*! Vector multiplication
    Vector& operator *= (const Vector& s)
    { *this = *this * s; return *this; }
  */

  //! Scalar multiplication
  template<typename U>
  Vector& operator *= (const U& a)
    { for (unsigned i=0; i<N; i++) x[i] *= a; return *this; }

  //! Scalar division
  template<typename U>
  Vector& operator /= (const U& a)
    { for (unsigned i=0; i<N; i++) x[i] /= a; return *this; }

  //! Equality
  bool operator == (const Vector& b) const
    { for(unsigned i=0; i<N; i++) if(x[i]!=b.x[i]) return false; return true; }

  //! Inequality
  bool operator != (const Vector& b) const
    { return ! operator==(b); }

  //! Vector addition
  template<typename U>
  const friend Vector operator + (Vector a, const Vector<N,U>& b)
    { a+=b; return a; }

  //! Vector subtraction
  template<typename U> 
  const friend Vector operator - (Vector a, const Vector<N,U>& b)
    { a-=b; return a; }

  //! Scalar multiplication
  const friend Vector operator * (Vector a, T c)
    { a*=c; return a; }

  //! Scalar multiplication
  const friend Vector operator * (T c, Vector a)
    { a*=c; return a; }

  //! Dot product
  const friend T operator * (const Vector& a, const Vector& b)
    { T r=0; for (unsigned i=0; i<N; i++) r += a[i]*b[i]; return r; }

  //! Scalar division
  const friend Vector operator / (Vector a, T c)
    { a/=c; return a; }

  //! Negation
  const friend Vector operator - (Vector s)
    { for (unsigned i=0; i<N; i++) s.x[i] = -s.x[i]; return s; }

  //! Access to elements
  T& operator [] (unsigned n)
    { return x[n]; }
  
  //! Alternative access to elements 
  const T operator [] (unsigned n) const
    { return x[n]; }

  //! Dimension of data
  unsigned size () const { return N; }

  //! Return the basis vector
  static const Vector basis (unsigned i)
    { Vector v;  v[i] = 1.0; return v; }

};

//! Cross product
template <typename T> 
const Vector<3,T> cross (const Vector<3,T>& a, const Vector<3,T>& b)
{
  Vector<3,T> result;
  unsigned j, k;
  for (unsigned i=0; i<3; i++) {
    j = (i+1)%3;  k = (i+2)%3;
    result[i] = a[j]*b[k] - a[k]*b[j];
  }

  return result;
}

template <typename T>
T normsq (const T& v)
{
  return v*v;
}

template <typename T>
T normsq (const std::complex<T>& v)
{
  return std::norm(v);
}

template <typename T, typename U=T>
T normsq (const class Estimate<T,U>& estimate)
{
  return normsq(estimate.get_value());
}

//! squared "length"/"norm"
template <unsigned N, typename T>
T normsq (const Vector< N, T>& v)
{
  T sum = normsq(v[0]);
  for (unsigned i=1; i < N; i++)
    sum += normsq(v[i]);
  return sum;
}

//! squared "length"/"norm"
template <unsigned N, template<typename T> class V, typename U >
U normsq (const Vector< N, V<U> >& v)
{
  U sum = normsq(v[0]);
  for (unsigned i=1; i < N; i++)
    sum += normsq(v[i]);
  return sum;
}

//! and the norm itself; note, won't work for complex
template <unsigned N, typename T>
T norm(const Vector<N, T> &v)
{
  return sqrt(normsq(v));
}


template<unsigned N, typename T>
Vector<N,T> real (const Vector< N, std::complex<T> >& input)
{
  Vector<N,T> result;
  for (unsigned i=0; i < N; i++)
    result[i] = input[i].real();
  return result;
}

template<unsigned N, typename T>
Vector<N,T> imag (const Vector< N, std::complex<T> >& input)
{
  Vector<N,T> result;
  for (unsigned i=0; i < N; i++)
    result[i] = input[i].imag();
  return result;
}

template<unsigned N, typename T>
Vector<N, std::complex<T> > conj (const Vector< N, std::complex<T> >& input)
{
  Vector< N, std::complex<T> > result;
  for (unsigned i=0; i < N; i++)
    result[i] = std::conj(input[i]);
  return result;
}

template<unsigned N, typename T>
bool myfinite (const Vector< N, T>& input)
{
  for (unsigned i=0; i < N; i++)
    if (!myfinite( input[i] ))
      return false;
  return true;
}

//! Maps the structure of Vector to other template methods
template<unsigned N, typename T> struct DatumTraits< Vector<N,T> >
{
  typedef T element_type;

  static inline unsigned ndim () { return N; }
  static inline T& element (Vector<N,T>& t, unsigned i) 
  { return t[i]; }
  static inline const T element (const Vector<N,T>& t, unsigned i)
  { return t[i]; }
};

//! Useful for quickly printing the components
template<unsigned N, typename T>
std::ostream& operator<< (std::ostream& ostr, const Vector<N,T>& v)
{
  ostr << "(" << v[0];
  for (unsigned i=1; i<N; i++)
    ostr << "," << v[i];
  return ostr << ")";
}

template<unsigned N, typename T>
std::istream& operator >> (std::istream& is, Vector<N,T>& v)
{
  char c;
  is >> c;

  if (c != '(') {
    is.setstate(std::istream::failbit);
    return is;
  }

  is >> v[0];

  for (unsigned i=1; i<N; i++)
  {
    is >> c >> v[i];
    if (c != ',') {
      is.setstate(std::istream::failbit);
      return is;
    }
  }

  is >> c;

  if (c != ')')
    is.setstate(std::istream::failbit);

  return is;
}

#endif  /* not __Vector_H defined */

