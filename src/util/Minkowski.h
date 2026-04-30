//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2014 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/Minkowski.h

#ifndef __epsic_util_Minkowski_H
#define __epsic_util_Minkowski_H

#include "Matrix.h"

//! Computes the Minkowski inner and outer products of two four-vectors
class Minkowski
{
public:

  //! Return the Minkowski inner product of two four-vectors
  template<typename T, typename U> typename PromoteTraits<T,U>::promote_type
  static inner (const Vector<4,T>& A, const Vector<4,U>& B)
  {
    typename PromoteTraits<T,U>::promote_type result = A[0] * B[0];
    for (unsigned i=1; i<4; i++)
      result -= A[i] * B[i];
    return result;
  }

  //! Return the Minkowski outer product of two four-vectors
  template<typename T, typename U> 
  Matrix<4,4,typename PromoteTraits<T,U>::promote_type>
  static outer (const Vector<4,T>& A, const Vector<4,U>& B)
  {
    Matrix<4,4,typename PromoteTraits<T,U>::promote_type> result = ::outer(A,B);
    typename PromoteTraits<T,U>::promote_type inv = 0.5 * inner(A,B);
    result[0][0] -= inv;
    for (unsigned i=1; i<4; i++)
      result[i][i] += inv;
    return result;
  }
};

#endif
