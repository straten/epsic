//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2004 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// psrchive/Util/units/Cloude.h

#ifndef __Cloude_H
#define __Cloude_H

#include "Pauli.h"
#include "Matrix.h"

//! Return the Hermitian target coherency matrix
/* This method computes the hermitian target coherency matrix, derived from
   the input Jones matrix as described by equations 4.6 and 4.9 of

   Cloude, S.R., 1986, "Group theory and polarization algebra",
   Optik, 75(1), 26-36
*/
template<typename T>
Matrix< 4,4, std::complex<T> > coherence (const Jones<T>& jones)
{
  // form the Hermitian biquaternion (Eq. 4.6)
  Quaternion<std::complex<T>,Hermitian> q = convert(jones);

  // convert to the equivalent complex vector and its Hermitian transpose
  Vector< 4, std::complex<T> > vect;
  Vector< 4, std::complex<T> > herm;

  for (unsigned i=0; i<4; i++) {
    vect[i] = q[i];
    herm[i] = conj(q[i]);
  }

  // Eq. 4.9
  return outer(vect,herm);
}


//! Return the Hermitian target coherency matrix derivative
template<typename T>
Matrix< 4,4, std::complex<T> > coherence (const Jones<T>& J,
					  const Jones<T>& Jgrad)
{
  // form the Hermitian biquaternion (Eq. 4.6)
  Quaternion<std::complex<T>,Hermitian> q = convert( J );
  Quaternion<std::complex<T>,Hermitian> qgrad = convert( Jgrad );

  // convert to the equivalent complex vector and its Hermitian transpose
  Vector< 4, std::complex<T> > v;
  Vector< 4, std::complex<T> > h;

  Vector< 4, std::complex<T> > vgrad;
  Vector< 4, std::complex<T> > hgrad;

  for (unsigned i=0; i<4; i++)
  {
    v[i] = q[i];
    h[i] = conj(q[i]);

    vgrad[i] = qgrad[i];
    hgrad[i] = conj(qgrad[i]);
  }

  return outer(v,hgrad) + outer(vgrad,h); 
}

//! Return the Jones matrix of the target coherency matrix left eigenvector
/* The left eigenvector, as returned by Jacobi, is the row vector
   given by the hermitian transpose of the right (column)
   eigenvector. */
template<typename T>
Jones<T> system (const Vector< 4,std::complex<T> >& left_eigen)
{
  Quaternion<std::complex<T>,Hermitian> q;
  
  for (unsigned idim=0; idim<4; idim++)
    q[idim] = conj( left_eigen[idim] );

  return convert (q);
}

//! Return the target entropy, given the target coherency matrix eigenvalues
/*! The eigenvalues are returned by Jacobi as a vector.  The target
  entropy is defined at the end of Section 4 of Cloude (1986). */
template<unsigned N, typename T>
  double entropy (const Vector<N,T>& eigenvalues)
{
  unsigned i;
  double total = 0;

  for (i=0; i<N; i++)
    total += eigenvalues[i];

  double entropy = 0;
  for (i=0; i<N; i++) {
    double Pi = eigenvalues[i]/total;
    entropy -= Pi * log(Pi)/log(double(N));
  }

  return entropy;
}

#endif
