//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2003 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/Pauli.h

#ifndef __Pauli_H
#define __Pauli_H

#include "Jones.h"
#include "Quaternion.h"
#include "Stokes.h"
#include "Basis.h"

#include <vector>
#include <limits>

namespace Pauli {

  //! The basis through which Stokes parameters are converted to Jones matrices
  Basis<double>& basis();

  //! Get the specified basis matrix
  Jones<double> matrix (unsigned i);
}


// convert Hermitian Quaternion to Jones matrix
template<typename T>
const Jones<T> convert (const Quaternion<std::complex<T>,Hermitian>& q)
{
  return Jones<T> (q.s0+q.s1, q.s2-ci(q.s3),
		   q.s2+ci(q.s3), q.s0-q.s1);
}

// convert Unitary Quaternion to Jones matrix
template<typename T>
const Jones<T> convert (const Quaternion<std::complex<T>,Unitary>& q)
{
  return Jones<T> (q.s0+ci(q.s1), q.s3+ci(q.s2),
		   -q.s3+ci(q.s2), q.s0-ci(q.s1));
}

// convert Hermitian Quaternion to Jones matrix
template<typename T>
const Jones<T> convert (const Quaternion<T,Hermitian>& q)
{
  return Jones<T> (std::complex<T>(q.s0+q.s1,0.0), std::complex<T>(q.s2,-q.s3),
		   std::complex<T>(q.s2,q.s3), std::complex<T>(q.s0-q.s1,0.0));
}

// convert Unitary Quaternion to Jones matrix
template<typename T>
const Jones<T> convert (const Quaternion<T,Unitary>& q)
{
  return Jones<T> (q.s0+ci(q.s1), q.s3+ci(q.s2),
		   -q.s3+ci(q.s2), q.s0-ci(q.s1));
}

// convert complex Stokes parameters to Jones matrix
template<typename T>
const Jones<T> convert (const Stokes< std::complex<T> >& stokes)
{
  Quaternion<std::complex<T>,Hermitian> q;

  q.set_scalar (stokes.get_scalar());
  q.set_vector (Pauli::basis().get_out(stokes.get_vector()));

  std::complex<T> half (0.5, 0.0);
  return convert (half*q);
}

// convert Stokes parameters to Jones matrix
template<typename T>
const Jones<T> convert (const Stokes<T>& stokes)
{
  Quaternion<T,Hermitian> q (stokes.get_scalar(),
                             Pauli::basis().get_out(stokes.get_vector()));
  return convert (T(0.5)*q);
}

// convert coherency vector to Jones matrix
template<typename T>
const Jones<T> convert (const std::vector<T>& c)
{
  if (c.size() != 4)
    throw std::runtime_error ("Jones<T> convert (std::vector<T>) vector.size() != 4");

  return Jones<T> (c[0], std::complex<T> (c[2], -c[3]),
                 std::complex<T> (c[2], c[3]), c[1]);
}

// convert Stokes parameters to natural basis
template<typename T>
const Quaternion<T,Hermitian> natural (const Stokes<T>& stokes)
{
  return Quaternion<T,Hermitian> 
    (stokes.get_scalar(), Pauli::basis().get_out(stokes.get_vector()));
}

template<typename T>
const Stokes<T> standard (const Quaternion<T,Hermitian>& q)
{ 
  return Stokes<T>( q.get_scalar(), Pauli::basis().get_in(q.get_vector()) );
}


// convert Jones matrix to Hermitian Biquaternion
template<typename T>
const Quaternion<std::complex<T>, Hermitian> convert (const Jones<T>& j)
{
  return Quaternion<std::complex<T>, Hermitian>
    ( T(0.5) * (j.j00 + j.j11),
      T(0.5) * (j.j00 - j.j11),
      T(0.5) * (j.j01 + j.j10),
      T(0.5) * ci (j.j01 - j.j10) );
}

// convert Jones matrix to Unitary Biquaternion
template<typename T>
const Quaternion<std::complex<T>, Unitary> unitary (const Jones<T>& j)
{
  return Quaternion<std::complex<T>, Unitary>
    ( T(0.5) *     (j.j00 + j.j11),
      T(-0.5) * ci (j.j00 - j.j11),
      T(-0.5) * ci (j.j01 + j.j10),
      T(0.5) *     (j.j01 - j.j10) );
}

// convert a Hermitian Jones matrix to Stokes parameters
template<typename T>
const Stokes<T> coherency (const Jones<T>& j)
{
  Quaternion<std::complex<T>,Hermitian> h = convert (j);
  return real_coherency (h);
}

// convert a Hermitian Jones matrix to Stokes parameters
template<typename T>
const Stokes< std::complex<T> > complex_coherency (const Jones<T>& j)
{
  Quaternion<std::complex<T>,Hermitian> h = convert (j);
  return coherency (h);
}

// convert a complex Hermitian Quaternion to Stokes parameters
template<typename T>
const Stokes<T> real_coherency (const Quaternion<std::complex<T>,Hermitian>& q)
{
  Quaternion<T,Hermitian> realpart (real(q));
  Quaternion<T,Hermitian> imaginary (imag(q));
  T nr = norm(realpart);
  T ni = norm(imaginary);
  if (ni > 1e8*std::numeric_limits<T>::epsilon() && ni*1e5 > nr)
  {
    const char* msg = "Stokes<T> coherency (Quaternion<complex<T>,Hermitian>) "
      "non-zero imaginary component";
    
    std::cerr << msg << 
      "\n   norm(imag(q))=" << ni << " norm=" << nr << std::endl;

    throw std::runtime_error (msg);
  }
  
  return coherency (realpart);
}

// convert a Hermitian Quaternion to Stokes parameters
template<typename T>
const Stokes<T> coherency (const Quaternion<T,Hermitian>& q)
{ 
  return Stokes<T>( T(2.0) * q.get_scalar(), 
		    T(2.0) * Pauli::basis().get_in(q.get_vector()) );
}

// transform the Stokes parameters by the given Jones matrix
template<typename T, typename U>
const Stokes<T> transform (const Stokes<T>& input, const Jones<U>& jones)
{
  return coherency (jones * convert(input) * herm(jones));
}

// transform the Stokes parameters by the given Jones matrix
template<typename T, typename U>
const Stokes< std::complex<T> >
transform (const Stokes< std::complex<T> >& input, const Jones<U>& jones)
{
  return complex_coherency (jones * convert(input) * herm(jones));
}

// transform the coherency matrix by the given Mueller matrix
template<typename T, typename U>
Jones<T> transform (const Matrix<4,4,U>& M, const Jones<T>& rho)
{
  Stokes< std::complex<T> > s = M * complex_coherency(rho);
  return convert(s);
}

// convert Jones matrix to Hermitian and Unitary Quaternion
template<typename T>
void polar (std::complex<T>& d, Quaternion<T, Hermitian>& h,
	    Quaternion<T, Unitary>& u, Jones<T> j)
{
  // make j unimodular
  d = sqrt (det(j));
  j /= d;

  // calculate the hermitian
  h = sqrt( real( convert( j*herm(j) ) ) );

  // remove the hermitian component from j
  j = inv(convert(h)) * j;

  // take the unitary component out of j
  u = real ( unitary (j) );
}


// multiply a Jones matrix by a Quaternion
template<typename T, typename U, QBasis B>
const Jones<T> operator * (const Jones<T>& j, const Quaternion<U,B>& q)
{
  return j * convert(q);
}

// multiply a Jones matrix by a Quaternion
template<typename T, typename U, QBasis B>
const Jones<T> operator * (const Quaternion<T,B>& q, const Jones<U>& j)
{
  return convert(q) * j;
}

// multiply Quaternions from different QBasis
template<typename T, typename U, QBasis A, QBasis B>
const Jones<T> operator * (const Quaternion<T,A>& q, const Quaternion<U,B>& u)
{
  return convert(q) * convert(u);
}

// return the Mueller matrix corresponding to the given Jones matrix
template<typename T>
Matrix<4,4,T> Mueller (const Jones<T>& J)
{
  Matrix<4,4,T> result;

  for (unsigned row=0; row < 4; row++) {
    Stokes<T> basis;
    basis[row] = 1.0;
    result[row] = coherency(herm(J) * convert(basis) * J);
  }

  return result;
}

// return the Mueller matrix derivative Jones matrix
template<typename T>
Matrix<4,4,T> Mueller (const Jones<T>& J, const Jones<T>& Jgrad)
{
  Matrix<4,4,T> result;

  for (unsigned row=0; row < 4; row++) {
    Stokes<T> basis;
    basis[row] = 1.0;
    Jones<T> rho = convert(basis);
    result[row] = coherency(herm(Jgrad) * rho * J + herm(J) * rho * Jgrad);
  }

  return result;
}

#endif
