//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2003 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/Matrix.h

#ifndef __Matrix_H
#define __Matrix_H

#include "Vector.h"
#include "complex_promote.h"
#include <stdexcept>

//! Matrix is a column vector of row vectors
template <unsigned Rows, unsigned Columns, typename T> 
class Matrix : public Vector< Rows, Vector<Columns,T> > 
{
public:

  //! Null constructor
  Matrix () { zero (); }

  //! Scalar constructor (identity times scalar)
  Matrix (T s) 
  {
    zero ();
    for (unsigned i=0; i<Rows; i++)
      this->x[i][i] = s;
  }

  //! Construct from another Vector of Vector<U> instance
  template<typename U> Matrix (const Vector< Rows, Vector<Columns,U> >& s)
  { operator=(s); }

  //! Set this instance equal to another Matrix<U> instance
  template<typename U> Matrix& operator = 
    (const Vector< Rows,Vector<Columns,U> >& s)
  {
    for (unsigned i=0; i<Rows; i++) 
      this->x[i] = s.x[i];
    return *this;
  }
  
  //! Negation
  const friend Matrix operator - (Matrix s)
  {
    for (unsigned i=0; i<Rows; i++) 
      for (unsigned j=0; j<Columns; j++)
	s.x[i][j] = -s.x[i][j];
    return s;
  }

  void zero ()
  {
    for (unsigned i=0; i<Rows; i++) 
      for (unsigned j=0; j<Columns; j++)
	this->x[i][j] = T(0.0);
  }
};


//! Vector multiplication
template<unsigned Rows, unsigned Columns, typename T, typename U>
const Vector<Rows,U> operator * (const Matrix<Rows,Columns,T>& m,
                                 const Vector<Columns,U>& b)
{
  Vector<Rows,U> r;
  for (unsigned i=0; i<Rows; i++)
    r[i] = Vector<Columns,U>(m[i]) * b;
  return r;
}

//! Vector transpose multiplication
template<unsigned Rows, unsigned Columns, typename T, typename U>
const Vector<Columns,U> operator * (const Vector<Rows,U>& b,
				    const Matrix<Rows,Columns,T>& m)
{
  Vector<Columns,U> r;
  for (unsigned j=0; j<Columns; j++)
    for (unsigned i=0; i<Rows; i++)
      r[j] += m[i][j] * b[i];
  return r;
}

//! Matrix multiplication
template<unsigned R1, unsigned C1R2, unsigned C2, typename T>
const Matrix<R1, C2, T>
operator * (const Matrix<R1, C1R2, T>& a, const Matrix<C1R2, C2, T>& b)
{ 
  Matrix<R1, C2, T> r; 
  for (unsigned i=0; i<R1; i++)
    for (unsigned j=0; j<C2; j++)
      for (unsigned k=0; k<C1R2; k++)
	r[i][j] += a[i][k]*b[k][j];
  return r;
}

template <unsigned Rows, unsigned C1, typename T, unsigned C2, typename U>
void GaussJordan (Matrix<Rows,C1,T>& a, Matrix<Rows,C2,U>& b)
{
  Vector<Rows, unsigned> indxc;
  Vector<Rows, unsigned> indxr;
  Vector<Rows, unsigned> ipiv;

  unsigned irow = 0;
  unsigned icol = 0;

  unsigned i, j, k;
  for (i=0; i<Rows; i++) {

    double big = 0.0;

    //cerr << "0" << endl;

    for (j=0; j<Rows; j++) {
      if (ipiv[j] != 1)
	for (k=0; k<Rows; k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  } 
	  else if (ipiv[k] > 1) 
	    throw std::runtime_error("GaussJordan (Matrix) Singular Matrix-1");
	}
    }

    //cerr << "1" << endl;

    ipiv[icol]++;
    
    if (irow != icol) {
      for (j=0; j<Rows; j++)
	std::swap (a[irow][j], a[icol][j]);
      for (j=0; j<C2; j++)
	std::swap (b[irow][j], b[icol][j]);
    }

    //cerr << "2" << endl;

    indxr[i]=irow;
    indxc[i]=icol;

    if (a[icol][icol] == 0.0)
      throw std::runtime_error ("GaussJordan (Matrix) Singular Matrix-2");

    //cerr << "3" << endl;

    T pivinv = 1.0/a[icol][icol];
    a[icol][icol]=1.0;

    for (j=0; j<Rows; j++)
      a[icol][j] *= pivinv;
    for (j=0; j<C2; j++)
      b[icol][j] *= pivinv;

    //cerr << "4" << endl;

    for (j=0; j<Rows; j++)
      if (j != icol) {
	T dum = a[j][icol];
	a[j][icol]=0.0;
	for (k=0; k<Rows; k++) 
	  a[j][k] -= a[icol][k]*dum;
	for (k=0; k<C2; k++)
	  b[j][k] -= b[icol][k]*dum;
      }

    //cerr << "5" << endl;
  }

  //cerr << "6" << endl;
    
  for (i=Rows; i>0; i--) {
    if (indxr[i-1] != indxc[i-1])
      for (j=0; j<Rows; j++)
	std::swap(a[j][indxr[i-1]],a[j][indxc[i-1]]);
  }
  
   //cerr << "7" << endl;
 
}

template <unsigned RC, typename T>
void matrix_identity (Matrix<RC, RC, T>& m)
{
  for (unsigned i=0; i<RC; i++)
    for (unsigned j=0; j<RC; j++)
      if (i==j)
	m[i][j] = 1;
      else
	m[i][j] = 0;
}

//! Returns the trace
template<unsigned Square, typename T>
T trace (const Matrix<Square,Square,T>& M)
{
  T result (0);
  for (unsigned i=0; i<Square; i++)
    result += M[i][i];
  return result;
}

template <unsigned RC, typename T>
const Matrix<RC, RC, T> inv (const Matrix<RC, RC, T>& m)
{
  Matrix<RC, RC, T> copy (m);
  
  Matrix<RC, RC, T> inverse;
  matrix_identity (inverse);

  GaussJordan (copy, inverse);

  return inverse;
}

template <unsigned Rows, unsigned Columns, typename T>
const Matrix<Columns, Rows, T> transpose (const Matrix<Rows, Columns,T>& m)
{
  Matrix< Columns, Rows,T> result;

  for (unsigned i=0; i<Rows; i++)
    for (unsigned j=0; j<Columns; j++)
      result[j][i] = m[i][j];

  return result;
}

template <unsigned Rows, unsigned Columns, typename T>
const Matrix< Columns, Rows,T> herm (const Matrix<Rows, Columns,T>& m)
{
  Matrix< Columns, Rows,T> result;

  for (unsigned i=0; i<Rows; i++)
    for (unsigned j=0; j<Columns; j++)
      result[j][i] = myconj( m[i][j] );

  return result;
}

//! Vector direct (outer) product 
template<unsigned Rows, unsigned Columns, typename T, typename U>
const Matrix<Rows,Columns,typename PromoteTraits<T,U>::promote_type> 
outer (const Vector<Rows,T>& a, const Vector<Columns,U>& b)
{
  Matrix<Rows, Columns, typename PromoteTraits<T,U>::promote_type> result;

  for (unsigned i=0; i<Rows; i++)
    for (unsigned j=0; j<Columns; j++)
      result[i][j] = a[i] * b[j];

  return result;
}

//! Matrix direct (Kronecker) product 
template<unsigned Ar, unsigned Ac, typename At,
	 unsigned Br, unsigned Bc, typename Bt>
const Matrix<Ar*Br,Ac*Bc,typename PromoteTraits<At,Bt>::promote_type> 
direct (const Matrix<Ar,Ac,At>& a, const Matrix<Br,Bc,Bt>& b)
{
  Matrix<Ar*Br,Ac*Bc, typename PromoteTraits<At,Bt>::promote_type> result;

  for (unsigned ar=0; ar<Ar; ar++)
    for (unsigned ac=0; ac<Ac; ac++)
      for (unsigned br=0; br<Br; br++)
	for (unsigned bc=0; bc<Bc; bc++)
	  result[ar*Br+br][ac*Bc+bc] = a[ar][ac] * b[br][bc];

  return result;
}

//! Partition a matrix into four regions
template<unsigned U, unsigned L, unsigned B, unsigned R, typename T>
void partition (const Matrix<U+B,L+R,T>& A,
		Matrix<U,L,T>& upper_left,
		Matrix<U,R,T>& upper_right,
		Matrix<B,L,T>& bottom_left,
		Matrix<B,R,T>& bottom_right)
{
  unsigned i, j;

  for (i=0; i<U; i++)
  {
    for (j=0; j<L; j++)
      upper_left[i][j] = A[i][j];
    for (j=0; j<R; j++)
      upper_right[i][j] = A[i][j+L];
  }

  for (i=0; i<B; i++)
  {
    for (j=0; j<L; j++)
      bottom_left[i][j] = A[i+U][j];
    for (j=0; j<R; j++)
      bottom_right[i][j] = A[i+U][j+L];
  }
}

//! Merge four matrices into one
template<unsigned U, unsigned L, unsigned B, unsigned R, typename T>
void compose (Matrix<U+B,L+R,T>& A,
	      const Matrix<U,L,T>& upper_left,
	      const Matrix<U,R,T>& upper_right,
	      const Matrix<B,L,T>& bottom_left,
	      const Matrix<B,R,T>& bottom_right)
{
  unsigned i, j;

  for (i=0; i<U; i++)
  {
    for (j=0; j<L; j++)
      A[i][j] = upper_left[i][j];
    for (j=0; j<R; j++)
      A[i][j+L] = upper_right[i][j];
  }

  for (i=0; i<B; i++)
  {
    for (j=0; j<L; j++)
      A[i+U][j] = bottom_left[i][j];
    for (j=0; j<R; j++)
      A[i+U][j+L] = bottom_right[i][j];
  }
}

//! Convenience interface for square, symmetric matrices
template<unsigned M, typename T>
void partition (const Matrix<M+1,M+1,T>& covariance,
		T& variance,
		Vector <M,T>& cov_vector,
		Matrix <M,M,T>& cov_matrix)
{
  Matrix<1,1,T> ul;
  Matrix<M,1,T> bl;
  Matrix<1,M,T> ur;

  partition (covariance, ul, ur, bl, cov_matrix);

  variance  = ul[0][0];
  cov_vector = ur[0];
}

//! Convenience interface for square, symmetric matrices
template<unsigned M, typename T>
void compose (Matrix<M+1,M+1,T>& covariance,
	      T variance,
	      const Vector <M,T>& cov_vector,
	      const Matrix <M,M,T>& cov_matrix)
{
  Matrix<1,1,T> ul;
  Matrix<M,1,T> bl;
  Matrix<1,M,T> ur;

  ul[0][0] = variance;

  for (unsigned i=0; i<M; i++)
    ur[0][i] = bl[i][0] = cov_vector[i];

  compose (covariance, ul, ur, bl, cov_matrix);
}

//! Return a 3-dimensional rotation about an axis
/*! After multiplying a Vector by the return value of this function,
    that Vector will be rotated about the axis, v, by an angle in
    radians in a direction that is consistent with the right-hand rule. */

template<typename T>
Matrix<3,3,T> rotation (const Vector<3,T>& v, double radians)
{
  double s = sin(radians);
  double c = cos(radians);
  double u = 1.0 - c;

  Matrix<3,3,T> result;

  result[0] = Vector<3,T>
    ( v[0]*v[0]*u + c    ,   v[1]*v[0]*u - v[2]*s,  v[2]*v[0]*u + v[1]*s );
  result[1] = Vector<3,T>
    ( v[0]*v[1]*u + v[2]*s,  v[1]*v[1]*u + c    ,   v[2]*v[1]*u - v[0]*s );
  result[2] = Vector<3,T>
    ( v[0]*v[2]*u - v[1]*s,  v[1]*v[2]*u + v[0]*s,  v[2]*v[2]*u + c );

  return result;
}

//! squared "length"/"norm"
template<unsigned R, unsigned C, typename T>
T normsq (const Matrix<R,C,T>& v)
{
  T sum = normsq(v[0]);
  for (unsigned i=1; i < R; i++)
    sum += normsq(v[i]);
  return sum;
}

//! Maps the structure of Matrix to other template methods
template<unsigned R, unsigned C, typename T> 
struct DatumTraits< Matrix<R,C,T> >
{
  typedef T element_type;

  static inline unsigned ndim () { return R*C; }
  static inline T& element (Matrix<R,C,T>& t, unsigned i) 
  { return t[i/C][i%C]; }
  static inline const T element (const Matrix<R,C,T>& t, unsigned i)
  { return t[i/C][i%C]; }
};


//! Useful for printing the components
template<unsigned R, unsigned C, typename T>
std::ostream& operator<< (std::ostream& ostr, const Matrix<R,C,T>& m)
{
  ostr << "[" << m[0];
  for (unsigned i=1; i<R; i++)
    ostr << ",\n " << m[i];
  return ostr << "]";
}

#endif  /* not __Matrix_H defined */

