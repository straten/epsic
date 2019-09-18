//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2004 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/Jacobi.h

#ifndef __Jacobi_H
#define __Jacobi_H

#include "Matrix.h"
#include "Quaternion.h"

template <typename T>
T norm (T x) { return fabs(x); }


//! Calculates the Jacobi rotation parameters for a real symmetric matrix
/*! This template function calculates the Jacobi transformation parameters
  that minimize roundoff error, as described in Section 11.1 of NR */
template <typename T, typename U>
void calculate_Jacobi (const U& p, const U& q, const T& pq,
                       T& s, T& tau, U& correction)
{
  T h = q - p;
  T g = 100.0 * fabs(pq);

  T t;
  
  if (fabs(h)+g == fabs(h))

    t=pq/h;

  else {

    // Equation 11.1.8
    T theta=0.5*h/(pq);

    // Equation 11.1.10
    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
    if (theta < 0.0) t = -t;

  }
  
  // Equation 11.1.11
  T c=1.0/sqrt(1+t*t);

  // Equation 11.1.12
  s=t*c;

  // Equation 11.1.18
  tau=s/(1.0+c);

  // last term of Equation 11.1.14 and 11.1.15
  correction = t*pq;

}
  
// returns the correction factor for the diagonal
template <typename T, typename U>
void calculate_Jacobi (const U& p, const U& q, const std::complex<T>& pq,
                       std::complex<T>& s, std::complex<T>& tau, U& correction)
{
  T Sq = 0.5 * (p-q);
  T Su = pq.real();
  T Sv = -pq.imag();

  Quaternion<T, Hermitian> herm (0, Sq, Su, Sv);
  Quaternion<T, Unitary> rotation = eigen (herm);

  T c = rotation.s0;
  s = std::complex<T> (-rotation.s3, -rotation.s2);

  // similar to Equation 11.1.18
  tau = conj(s)/(1.0+c);

  correction = 2.0 * ( c * ( s*conj(pq) ).real() + Sq * ( s*conj(s) ).real() );
}


//! Perform the Jacobi rotation
/*! This template function performs the Jacobi transformation that
  minimizes roundoff error, as described in Section 11.1 of NR */
template <unsigned RC, typename T>
inline void rotate_Jacobi (Matrix<RC,RC,T>& x, T s, T tau,
			   unsigned i, unsigned j, unsigned k, unsigned l)
{
  T g=x[i][j];
  T h=x[k][l];

  // Equation 11.1.16
  x[i][j] -= myconj(s)*(h+g*myconj(tau));
  // Equation 11.1.17
  x[k][l] += s*(g-h*tau);
}

// returns the correction factor for the diagonal
template <unsigned RC, typename T, typename U>
U JacobiRotation ( unsigned ip, unsigned iq, 
		   Matrix<RC,RC,T>& a, Matrix<RC,RC,T>& v, Vector<RC,U>& d )
{
  T s, tau;
  U correction;

  //cerr << "calculate_Jacobi ip=" << ip << " iq=" << iq << endl;
  calculate_Jacobi (d[ip], d[iq], a[ip][iq], s, tau, correction);

  d[ip] -= correction;
  d[iq] += correction;
  
  unsigned j;

  //cerr << "in: pp=" << a[ip][ip] << " qq=" << a[iq][iq] << endl;

#if 0

  for (j=0; j<ip; j++)
    rotate_Jacobi(a,s,tau,j,ip,j,iq);
  for (j=ip+1; j<iq; j++)
    rotate_Jacobi(a,s,tau,ip,j,j,iq);
  for (j=iq+1; j<RC; j++)
    rotate_Jacobi(a,s,tau,ip,j,iq,j);

  a[ip][iq]=0.0;

#else

  // rotate the columns of the input matrix
  for (j=0; j<RC; j++)
    rotate_Jacobi(a,s,tau,j,ip,j,iq);

  s = myconj(s);
  tau = myconj(tau);

  // rotate the rows of the input matrix
  for (j=0; j<RC; j++)
    rotate_Jacobi(a,s,tau,ip,j,iq,j);

  //cerr << "out: pp=" << a[ip][ip] << " qq=" << a[iq][iq] << endl;
  //cerr << "d: pp=" << d[ip] << " qq=" << d[iq] << endl;

  a[iq][ip]=a[ip][iq]=0.0;

#endif
  
  // rotate the rows of the eigenvector matrix
  for (j=0; j<RC; j++)
    rotate_Jacobi (v,s,tau,ip,j,iq,j);

  return correction;
}



template <typename T, typename U, unsigned RC>
void Jacobi (Matrix<RC,RC,T>& a, Matrix<RC,RC,T>& evec, Vector<RC,U>& eval)
{
  Vector<RC,U> b, z;

  // copy the diagonal elements into d and b
  unsigned ip, iq;
  for (ip=0; ip<RC; ip++)
    b[ip] = eval[ip] = myreal( a[ip][ip] );

  // start with the identity matrix
  matrix_identity (evec);

  for (unsigned iter=0; iter < 50; iter++) {

    // sum off diagonal elements
    U sum = 0.0;

    for (ip=0; ip < RC; ip++)
      for (iq=ip+1; iq < RC; iq++)
	sum += norm(a[ip][iq]);

    if (sum == 0.0)
      return;

    U thresh = 0;

    if (iter < 4)
      thresh=0.2*sum/(RC*RC);

    for (ip=0; ip < RC; ip++)
      for (iq=ip+1; iq < RC; iq++) {

	// after four iterations, skip the rotation if the off-diagonal
	// element is small

	U g = 100.0 * norm(a[ip][iq]);

	if ( iter > 4 
	     && norm(eval[ip])+g == norm(eval[ip])
	     && norm(eval[iq])+g == norm(eval[iq]) )

	  a[ip][iq]=a[iq][ip]=0.0;


	else if (norm(a[ip][iq]) > thresh) {

	  U correction = JacobiRotation ( ip, iq, a, evec, eval );

	  z[ip] -= correction;
	  z[iq] += correction;

	}
      }

    for (ip=0; ip < RC; ip++) {
      b[ip] += z[ip];
      eval[ip]=b[ip];
      z[ip]=0.0;
    }


  } // for each iteration
  
}


#endif
