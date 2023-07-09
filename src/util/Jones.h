//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2003 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/Jones.h

#ifndef __Jones_H
#define __Jones_H

#include "Matrix.h"
#include "Traits.h"
#include "complex_promote.h"
#include "complex_math.h"

//! Jones matrices are 2x2 matrices with complex elements
template<typename T> class Jones {
  
public:
  std::complex<T> j00,j01,j10,j11;

  //! Default constructor
  Jones (T scalar = 0.0)
    : j00(scalar), j01(0.0), j10(0.0), j11(scalar) { }

  //! Construct from a castable type
  template<typename U> explicit Jones (const U& scalar)
    : j00(scalar), j01(0.0), j10(0.0), j11(scalar) { }

  //! Construct from std::complex<T>
  Jones (std::complex<T> j00_, std::complex<T> j01_,
	 std::complex<T> j10_, std::complex<T> j11_)
    : j00(j00_), j01(j01_), j10(j10_), j11(j11_) {  }

  //! Construct from another Jones<T> matrix
  Jones (const Jones& s)
    : j00(s.j00), j01(s.j01), j10(s.j10), j11(s.j11) {  }

  //! Construct from another Jones<U> matrix
  template<typename U> Jones (const Jones<U>& s)
    { operator=(s); }

  //! Construct from a Matrix
  template<typename U> Jones (const Matrix< 2, 2, std::complex<U> >& M)
    : j00(M[0][0]), j01(M[0][1]), j10(M[1][0]),j11(M[1][1]) {  }

  //! Set this instance equal to another Jones<T> instance
  Jones& operator = (const Jones& s)
    { j00=s.j00; j01=s.j01; j10=s.j10; j11=s.j11; return *this; }

  //! Set this instance equal to a scalar
  Jones& operator = (T scalar)
    { j00=scalar; j01=0; j10=0; j11=scalar; return *this; }

  //! Set this instance equal to a complex scalar
  Jones& operator = (std::complex<T> scalar)
    { j00=scalar; j01=0; j10=0; j11=scalar; return *this; }

  //! Set this instance equal to another Jones<U> instance
  template<typename U> Jones& operator = (const Jones<U>& s)
    { j00=std::complex<T>(s.j00.real(), s.j00.imag()); 
      j01=std::complex<T>(s.j01.real(), s.j01.imag());
      j10=std::complex<T>(s.j10.real(), s.j10.imag()); 
      j11=std::complex<T>(s.j11.real(), s.j11.imag()); return *this; }

  typedef Matrix< 2, 2, std::complex<T> > equiv;

  //! Cast to Matrix
  operator equiv () const
  { equiv M; M[0][0]=j00; M[0][1]=j01; M[1][0]=j10; M[1][1]=j11; return M; }

  //! Add another Jones<T> instance to this one
  Jones& operator += (const Jones& s)
    { j00+=s.j00; j01+=s.j01; j10+=s.j10; j11+=s.j11; return *this; }

  //! Subract another Jones<T> instance from this one
  Jones& operator -= (const Jones& s)
    { j00-=s.j00; j01-=s.j01; j10-=s.j10; j11-=s.j11; return *this; }

  //! Multiply another Jones<T> instance into this one (this=this*j)
  Jones& operator *= (const Jones& j);

  //! Multiply this instance by std::complex<U>
  template<typename U> Jones& operator *= (const std::complex<U>& au)
    { std::complex<T>a(au); j00*=a; j01*=a; j10*=a; j11*=a; return *this; }

  //! Divide this instance by std::complex<U>
  template<typename U> Jones& operator /= (const std::complex<U>& au)
    { std::complex<T>a(T(1.0),T(0.0));
      a/=au; j00*=a; j01*=a; j10*=a; j11*=a; return *this; }

  //! Multiply this instance by T
  Jones& operator *= (T a)
    { j00*=a; j01*=a; j10*=a; j11*=a; return *this; }

  //! Divide this instance by T
  Jones& operator /= (T a)
    { T d=1.0/a; j00*=d; j01*=d; j10*=d; j11*=d; return *this; }

  //! Equality
  bool operator == (const Jones& b) const
  { return 
      j00 == b.j00  &&  j01 == b.j01 && 
      j10 == b.j10  &&  j11 == b.j11;
  }

  //! Inequality
  bool operator != (const Jones& b) const
  { return ! Jones::operator==(b); }

  //! Binary multiplication of Jones<T> and std::complex<U>
  template<typename U> 
  const friend Jones operator * (Jones a, std::complex<U> c)
    { a*=c; return a; }

  //! Binary multiplication of std::complex<U> and Jones<T>
  template<typename U>
  const friend Jones operator * (std::complex<U> c, Jones a)
    { a*=c; return a; }

  //! Binary multiplication of Jones<T> and T
  const friend Jones operator * (Jones a, T c)
    { a*=c; return a; }

  //! Binary multiplication of T and Jones<T>
  const friend Jones operator * (T c, Jones a)
    { a*=c; return a; }

  //! Binary division of Jones by any type
  template<typename U> const friend Jones operator / (Jones a, U c)
    { a/=c; return a; }

  //! Negation operator returns negative of instance
  const friend Jones operator - (Jones s)
    { s.j00=-s.j00; s.j01=-s.j01; s.j10=-s.j10; s.j11=-s.j11; return s; }

  //! Returns reference to the value of the matrix at j(ir,ic)
  std::complex<T>& operator () (unsigned ir, unsigned ic)
  { std::complex<T>* val = &j00; return val[ir*2+ic]; }
  
  //! Returns const reference to the value of the matrix at j(ir,ic)
  const std::complex<T>& operator () (unsigned ir, unsigned ic) const
    { const std::complex<T>* val = &j00; return val[ir*2+ic]; }

  //! Alternative access to elements
  std::complex<T>&  operator [] (unsigned n)
  { std::complex<T>* val = &j00; return val[n]; }

  //! Alternative access to elements 
  const std::complex<T>& operator [] (unsigned n) const
  { const std::complex<T>* val = &j00; return val[n]; }

  //! Return true if the off-diagonal elements are zero
  const bool is_diagonal () const { return (j01==std::complex<T>(0.0)) 
    && (j10==std::complex<T>(0.0)); }

  //! The identity matrix
  static const Jones& identity();

  //! Dimension of data
  unsigned size () const { return 4; }

  //! Degree of polarization
  T p () const
  { 
    T tr = trace(*this).real();
    T d = det(*this).real();
    return sqrt( 1 - 4*d/(tr*tr) );
  }

};

//! Binary addition
template<typename T, typename U>
const Jones<typename PromoteTraits<T,U>::promote_type>
operator + (const Jones<T>& a, const Jones<U>& b)
{ Jones<typename PromoteTraits<T,U>::promote_type> ret(a); ret+=b; return ret; }

//! Binary subtraction
template<typename T, typename U>
const Jones<typename PromoteTraits<T,U>::promote_type> 
operator - (const Jones<T>& a, const Jones<U>& b)
{ Jones<typename PromoteTraits<T,U>::promote_type> ret(a); ret-=b; return ret; }

//! Binary multiplication of two Jones matrices
template<typename T, typename U>
const Jones<typename PromoteTraits<T,U>::promote_type> 
operator * (const Jones<T>& a, const Jones<U>& b)
{ Jones<typename PromoteTraits<T,U>::promote_type> ret(a); ret*=b; return ret; }


//! The identity matrix
template<typename T>
const Jones<T>& Jones<T>::identity ()
{
  static Jones<T> I (1,0,
 	             0,1);
  return I;
}

//! Enable the Jones class to be passed to certain template functions
template<typename T> struct DatumTraits< Jones<T> >
{
  ElementTraits< std::complex<T> > element_traits;
  typedef std::complex<T> element_type;

  static inline unsigned ndim () { return 4; }
  static inline std::complex<T>& element (Jones<T>& t, unsigned i) 
  { return t[i]; }
  static inline const std::complex<T>& element (const Jones<T>& t, unsigned i)
  { return t[i]; }
};

//! Multiply another Jones<T> instance into this one (this=this*j)
template<typename T>
Jones<T>& Jones<T>::operator *= (const Jones<T>& j)
{
  std::complex<T> temp (j00 * j.j00 + j01 * j.j10);
  j01  = j00 * j.j01 + j01 * j.j11; j00=temp;
  temp = j10 * j.j00 + j11 * j.j10;
  j11  = j10 * j.j01 + j11 * j.j11; j10=temp;
  return *this; 
}

//! Returns the inverse
template<typename T>
Jones<T> inv (const Jones<T>& j)
{
  std::complex<T> d(1.0,0.0); d/=det(j);
  return Jones<T>(d*j.j11, -d*j.j01,
		  -d*j.j10, d*j.j00);
}

//! Returns the complex conjugate
template<typename T>
Jones<T> conj (const Jones<T>& j)
{
  return Jones<T>(std::conj(j.j00), std::conj(j.j01),
		  std::conj(j.j10), std::conj(j.j11));
}

//! Returns the Hermitian transpose (transpose of complex conjugate)
template<typename T> 
Jones<T> herm (const Jones<T>& j)
{
  return Jones<T>(std::conj(j.j00), std::conj(j.j10),
		  std::conj(j.j01), std::conj(j.j11));
}

//! Returns the determinant
template<typename T>
std::complex<T> det (const Jones<T>& j) { return j.j00*j.j11 - j.j01*j.j10; }

//! Returns the trace
template<typename T>
std::complex<T> trace (const Jones<T>& j) { return j.j00 + j.j11; }

//! Returns the variance (square of the Frobenius norm)
template<typename T>
T norm (const Jones<T>& j)
{ return
    norm(j.j00) + norm(j.j01) + 
    norm(j.j10) + norm(j.j11);
}

template<typename T>
bool myfinite (const Jones<T>& j)
{ return myfinite(j.j00) && myfinite(j.j01) && myfinite(j.j10) && myfinite(j.j11); }

template<typename T>
T fabs (const Jones<T>& j)
{ 
  return sqrt (norm(j));
}

//! Useful for quickly printing the values of matrix elements
template<typename T>
std::ostream& operator<< (std::ostream& ostr, const Jones<T>& j)
{
  return ostr << "[" << j.j00 << j.j01 << j.j10 << j.j11 << "]";
}

#endif  /* not __Jones_H defined */

