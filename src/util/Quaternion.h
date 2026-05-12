//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2003-2009 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/Quaternion.h

#ifndef __epsic_util_Quaternion_H
#define __epsic_util_Quaternion_H

#include "Traits.h"
#include "complex_promote.h"
#include "Vector.h"

//! The type of quaternion determines if the vector part is real or imaginary
/*!
  Although the product of two Hermitian matrices is not necessarily Hermitian
  (and therefore is not a group operation), it is still useful to represent 
  Hermitian matrices in scalar+vector form.  The Quaternion class is used;
  however, calling something a Hermitian quaternion is conceptually misleading.
*/
enum QType { Hermitian, Unitary };

//! Represents either a left-handed quaternion or a biquaternion
/*!
  The Quaternion template class was designed to represent either biquaternions
  (a quaternion with complex-valued coefficients) or left-handed quaternions 
  using the %Pauli spin matrices as the basis (see below).
  
  Left-handed quaternions can have
  either purely imaginary (QType = Unitary) or purely real (QType = Hermitian)
  vector components. However, as the product of two Hermitian matrices is not
  Hermitian (unless the matrices commute), it is not possible to
  multiply two real quaternions in the %Pauli basis.

  There are four ways to create a null-constructed instance
  of a Quaternion object:

    \code{.cpp}
    // Left-handed quaternion (unitary basis)
    Quaternion<float> q;

    // Representation of Hermitian matrix (non-multiplicative)
    Quaternion<float, Hermitian> h;

    // Biquaternion in unitary basis
    Quaternion<complex<float>, Unitary> bq;

    // Biquaternion in Hermitian basis
    Quaternion<complex<float>, Hermitian> bh;
    \endcode

  Note that the floating-point precision is specified in the first
  template argument and the type, Unitary or Hermitian,
  is specified in the second template argument. As the left-handed
  quaternion is most closely related to Hamilton's quaternion, the
  Unitary basis is used by default.

  A biquaternion may be represented by the linear combination,
  \f[
    \bm{B}=b_0\pauli{0}+b_1\pauli{1}+b_2\pauli{2}+b_3\pauli{3},\hspace{5mm}b_i\in\C
  \f]
  where \f$\pauli{0}\f$ is the identity matrix and \f$\pauli{1-3}\f$
  are the Pauli spin matrices. The Pauli matrices have the following properties:
  \f[
    \pauli{i}^2 = \bm{I}, \hspace{5mm}
    \pauli{i}\pauli{j} = -\pauli{j}\pauli{i} = \Ci\pauli{k}
  \f]
  where \f$\{i,j,k\}\f$ is chosen from cyclic permutations of \f$\{1,2,3\}\f$.

  A biquaternion may be written as a scalar plus a three-vector:
  \f[
    \bm{B}=[b+\bm{b}] = b\pauli{0} + \bm{b\cdot\sigma}, \hspace{5mm} b_i\in\C
  \f]
  where \f$\bm{b}=(b_1,b_2,b_3)\f$ and \f$\bm{\sigma}\f$ is a three-vector whose
  components are the Pauli spin matrices.  Using this notation, the left-handed
  quaternion sub-group may be represented as
  \f[
    \bm{Q}=[q+i\bm{q}]= q\pauli{0} + \Ci\bm{q\cdot\sigma}, \hspace{5mm} q_i\in\R.
  \f]

  Consider the product of two biquaternions, \f$\bm{A}=[a_0+\bm{a}]\f$
  and \f$\bm{B}=[b_0+\bm{b}]\f$, 
  \f[
    \begin{array}{ccl}
    \bm{A}\bm{B} & = &
    (a_0\pauli{0} + a_1\pauli{1} + a_2\pauli{2} + a_3\pauli{3})
    (b_0\pauli{0} + b_1\pauli{1} + b_2\pauli{2} + b_3\pauli{3}) \\
    & = & a_0b_0\pauli{0}\pauli{0} + a_0b_1\pauli{0}\pauli{1}
            + a_0b_2\pauli{0}\pauli{2} + a_0b_3\pauli{0}\pauli{3} \\
    & + & a_1b_0\pauli{1}\pauli{0} + a_1b_1\pauli{1}\pauli{1}
            + a_1b_2\pauli{1}\pauli{2} + a_1b_3\pauli{1}\pauli{3} \\
    & + & a_2b_0\pauli{2}\pauli{0} + a_2b_1\pauli{2}\pauli{1}
            + a_2b_2\pauli{2}\pauli{2} + a_2b_3\pauli{2}\pauli{3} \\
    & + & a_3b_0\pauli{3}\pauli{0} + a_3b_1\pauli{3}\pauli{1}
            + a_3b_2\pauli{3}\pauli{2} + a_3b_3\pauli{3}\pauli{3}
    \end{array}
  \f]
  which reduces to
  \f[
    \begin{array}{ccl}
    \bm{A}\bm{B}
    & = & (a_0b_0 + a_1b_1 + a_2b_2 + a_3b_3) \pauli{0} \\
    & + & (a_0b_1 + a_1b_0 + \Ci a_2b_3 - \Ci a_3b_2) \pauli{1} \\
    & + & (a_0b_2 - \Ci a_1b_3 + a_2b_0 + \Ci a_3b_1) \pauli{2} \\
    & + & (a_0b_3 + \Ci a_1b_2 - \Ci a_2b_1 + a_3b_0) \pauli{3}
    \end{array}
  \f]

  This equation is used to derive the multiplication rule for 
  left-handed quaternions.  Assuming that
  \f$\bm{C}=[c_0+\Ci\bm{c}]\f$ and \f$\bm{D}=[d_0+\Ci\bm{d}]\f$,
  the left-handed quaternion product, \f$\bm{X}=[x_0+\Ci\bm{x}]=\bm{CD}\f$, 
  is given by
  \f[
    \begin{array}{rl}
    x_0 & = c_0d_0 - c_1d_1 - c_2d_2 - c_3d_3 \\
    x_1 & = c_0d_1 + c_1d_0 - c_2d_3 + c_3d_2 \\
    x_2 & = c_0d_2 + c_1d_3 + c_2d_0 - c_3d_1 \\
    x_3 & = c_0d_3 - c_1d_2 + c_2d_1 + c_3d_0
    \end{array}
  \f]

  Notice that this result differs from the product of two quaternions as
  originally defined by Hamilton.  This is because the left-handed
  quaternion basis is given by the Pauli spin matrices multiplied by the
  imaginary number, \f$\Ci=\sqrt{-1}\f$, which satisfy
  \f[
    \Ci\pauli{i} \Ci\pauli{j} = -\Ci\pauli{k}
  \f]

*/
template<typename T, QType B = Unitary> 
class Quaternion {
  
public:
  T s0,s1,s2,s3;

  //! Default constructor
  Quaternion (const T& a=0.0, const T& b=0.0, const T& c=0.0, const T& d=0.0)
    : s0(a), s1(b), s2(c), s3(d) { }

  //! Construct from a scalar and vector
  template<typename U> Quaternion (T s, const Vector<3, U>& v)
    : s0(s), s1(v[0]), s2(v[1]), s3(v[2]) { }

  //! Construct from another Quaternion<U> instance
  template<typename U> Quaternion (const Quaternion<U, B>& s)
    : s0(s.s0), s1(s.s1), s2(s.s2), s3(s.s3) { }

  //! Set this instance equal to another Quaternion<U> instance
  template<typename U>
  const Quaternion& operator = (const Quaternion<U,B>& s)
    { s0=T(s.s0); s1=T(s.s1); s2=T(s.s2); s3=T(s.s3); return *this; }

  //! Scalar addition
  const Quaternion& operator += (const T& s)
    { s0+=s; return *this; }

  //! Scalar subtraction
  const Quaternion& operator -= (const T& s)
    { s0-=s; return *this; }

  //! Scalar multiplication
  const Quaternion& operator *= (const T& a)
    { s0*=a; s1*=a; s2*=a; s3*=a; return *this; }

  //! Scalar division
  const Quaternion& operator /= (const T& a)
    { T d(1.0); d/=a; s0*=d; s1*=d; s2*=d; s3*=d; return *this; }

  //! Quaternion addition
  template<class U>
  const Quaternion& operator += (const Quaternion<U,B>& s)
    { s0+=s.s0; s1+=s.s1; s2+=s.s2; s3+=s.s3; return *this; }

  //! Quaternion subtraction
  template<class U>
  const Quaternion& operator -= (const Quaternion<U,B>& s)
    { s0-=s.s0; s1-=s.s1; s2-=s.s2; s3-=s.s3; return *this; }

  //! Quaternion multiplication
  template<class U>
  const Quaternion& operator *= (const Quaternion<U,B>& s)
    { *this = *this * s; return *this; }

  //! Equality
  bool operator == (const Quaternion& b) const
    { return  s0==b.s0 && s1==b.s1 && s2==b.s2 && s3==b.s3; }

  //! Equality
  bool operator == (T scalar) const
    { return  s0==scalar && s1==0 && s2==0 && s3==0; }

  //! Inequality
  bool operator != (const Quaternion& b) const
    { return ! operator==(b); }

  //! Negation
  friend Quaternion operator - (Quaternion s)
    { s.s0=-s.s0; s.s1=-s.s1; s.s2=-s.s2; s.s3=-s.s3; return s; }

  //! Alternative access to elements
  T& operator [] (int n)
    { T* val = &s0; return val[n]; }
  
  //! Alternative access to elements 
  T operator [] (int n) const
    { return *(&s0+n); }

  //! Access to scalar component
  T get_scalar () const { return s0; }

  //! Set the scalar component
  void set_scalar (T s) { s0 = s; }

  //! Access to vector component
  Vector<3,T> get_vector () const 
    { Vector<3,T> ret; ret[0]=s1; ret[1]=s2; ret[2]=s3; return ret; }

  //! Set the vector component
  template<typename U>
  void set_vector (const Vector<3,U>& v) { s1=v[0]; s2=v[1]; s3=v[2]; }

  //! Identity
  static const Quaternion& identity();

  //! Dimension of data
  unsigned size () const { return 4; }

};

//! Quaternion addition
template<typename T, QType B, typename U>
Quaternion<typename PromoteTraits<T,U>::promote_type,B>
operator + (const Quaternion<T,B>& a, const Quaternion<U,B>& b)
{
  Quaternion<typename PromoteTraits<T,U>::promote_type,B> ret(a);
  ret+=b;
  return ret;
}

//! Quaternion subtraction
template<typename T, QType B, typename U>
Quaternion<typename PromoteTraits<T,U>::promote_type,B>
operator - (const Quaternion<T,B>& a, const Quaternion<U,B>& b)
{
  Quaternion<typename PromoteTraits<T,U>::promote_type,B> ret(a);
  ret-=b;
  return ret;
}

//! Scalar multiplication
/*! The return type should use PromoteTraits, but the compiler won't */
template<typename T, QType B, typename U>
Quaternion<T,B>
operator * (const Quaternion<T,B>& a, const U& c)
{ 
  Quaternion<typename PromoteTraits<T,U>::promote_type,B> ret(a);
  ret*=c;
  return ret;
}

//! Scalar multiplication
/*! The return type should use PromoteTraits, but the compiler won't */
template<typename T, QType B, typename U>
Quaternion<T,B>
operator * (const U& c, const Quaternion<T,B>& a)
{
  Quaternion<typename PromoteTraits<T,U>::promote_type,B> ret(a);
  ret*=c;
  return ret;
}

//! Scalar division
template<typename T, QType B, typename U>
Quaternion<typename PromoteTraits<T,U>::promote_type,B>
operator / (const Quaternion<T,B>& a, const U& c)
{
  Quaternion<typename PromoteTraits<T,U>::promote_type,B> ret(a);
  ret/=c;
  return ret;
}

//! The identity Quaternion
template<typename T, QType B>
const Quaternion<T,B>& Quaternion<T,B>::identity ()
{
  static Quaternion<T,B> I (1,0,0,0);
  return I;
}

//! Enable the Quaternion class to be passed to certain template functions
template<typename T, QType B> struct DatumTraits< Quaternion<T,B> >
{
  ElementTraits<T> element_traits;
  static inline unsigned ndim () { return 4; }
  static inline T& element (Quaternion<T,B>& t, unsigned idim) 
  { return t[idim]; }
  static inline const T& element (const Quaternion<T,B>& t, unsigned idim)
  { return t[idim]; }
};

//! Multiplication of two Biquaternions in the Hermitian basis
template<typename T, typename U>
Quaternion<std::complex<typename PromoteTraits<T,U>::promote_type>, Hermitian>
operator * (const Quaternion<std::complex<T>,Hermitian>& a,
	    const Quaternion<std::complex<U>,Hermitian>& b)
{
  typedef std::complex<typename PromoteTraits<T,U>::promote_type> R;
  return Quaternion<R,Hermitian>
    ( a.s0*b.s0 + a.s1*b.s1 + a.s2*b.s2 + a.s3*b.s3 ,
      a.s0*b.s1 + a.s1*b.s0 + ci(a.s2*b.s3) - ci(a.s3*b.s2) ,
      a.s0*b.s2 - ci(a.s1*b.s3) + a.s2*b.s0 + ci(a.s3*b.s1) ,
      a.s0*b.s3 + ci(a.s1*b.s2) - ci(a.s2*b.s1) + a.s3*b.s0 );
}

//! Multiplication of two Quaternions in the Unitary basis
template<typename T, typename U>
Quaternion<typename PromoteTraits<T,U>::promote_type, Unitary>
operator * (const Quaternion<T,Unitary>& a, const Quaternion<U,Unitary>& b)
{
  return Quaternion<typename PromoteTraits<T,U>::promote_type, Unitary>
    (a.s0*b.s0 - a.s1*b.s1 - a.s2*b.s2 - a.s3*b.s3,
     a.s0*b.s1 + a.s1*b.s0 - a.s2*b.s3 + a.s3*b.s2,
     a.s0*b.s2 + a.s1*b.s3 + a.s2*b.s0 - a.s3*b.s1,
     a.s0*b.s3 - a.s1*b.s2 + a.s2*b.s1 + a.s3*b.s0);
}


//! Returns the real component of a Biquaternion
template<typename T, QType B>
Quaternion<T,B> real (const Quaternion<std::complex<T>,B>& j)
{
  return Quaternion<T,B>
    (j.s0.real(), j.s1.real(), j.s2.real(), j.s3.real());
}

//! Returns the imag component of a Biquaternion
template<typename T, QType B>
Quaternion<T,B> imag (const Quaternion<std::complex<T>,B>& j)
{
  return Quaternion<T,B>
    (j.s0.imag(), j.s1.imag(), j.s2.imag(), j.s3.imag());
}

//! Returns the complex conjugate of a Hermitian Biquaternion
template<typename T>
Quaternion<T,Hermitian> conj (const Quaternion<T,Hermitian>& j)
{
  return Quaternion<T,Hermitian>
    (myconj(j.s0), myconj(j.s1), myconj(j.s2), -myconj(j.s3));
}

//! Returns the complex conjugate of a Unitary Quaternion
template<typename T>
Quaternion<T,Unitary> conj (const Quaternion<T,Unitary>& j)
{ 
  return Quaternion<T,Unitary>
    (myconj(j.s0), -myconj(j.s1), -myconj(j.s2), myconj(j.s3));
}


//! Returns the Hermitian transpose of a Hermitian Quaternion
template<typename T>
Quaternion<T, Hermitian> herm (const Quaternion<T,Hermitian>& j)
{
  return Quaternion<T,Hermitian>
    (myconj(j.s0), myconj(j.s1), myconj(j.s2), myconj(j.s3));
}


//! Returns the Hermitian transpose of a Unitary Quaternion
template<typename T>
Quaternion<T, Unitary> herm (const Quaternion<T,Unitary>& j)
{
  return Quaternion<T,Unitary>
    (myconj(j.s0), -myconj(j.s1), -myconj(j.s2), -myconj(j.s3));
}


//! Returns the inverse of Quaternion, j
template<typename T, QType B>
Quaternion<T, B> inv (const Quaternion<T,B>& j) 
{
  T d (-1.0); d/=det(j);
  return Quaternion<T,B> (-d*j.s0, d*j.s1, d*j.s2, d*j.s3);
}


//! Returns the determinant of a Hermitian Quaternion
template<typename T>
T det (const Quaternion<T,Hermitian>& j)
{ return j.s0*j.s0 - j.s1*j.s1 - j.s2*j.s2 - j.s3*j.s3; }

//! Returns the determinant of a Unitary Quaternion
template<typename T>
T det (const Quaternion<T,Unitary>& j)
{ return j.s0*j.s0 + j.s1*j.s1 + j.s2*j.s2 + j.s3*j.s3; }

//! Returns the trace of Quaternion, j
template<typename T, QType B>
T trace (const Quaternion<T,B>& j)
{ return 2.0 * j.s0; }

//! Returns the square of the Frobenius norm of a Biquaternion
template<typename T, QType B>
T norm (const Quaternion<std::complex<T>,B>& j)
{ return 2.0 * (norm(j.s0) + norm(j.s1) + norm(j.s2) + norm(j.s3)); }

//! Returns the square of the Frobenius norm of a quaternion
template<typename T, QType B>
T norm (const Quaternion<T,B>& j)
{ return 2.0 * (j.s0*j.s0 + j.s1*j.s1 + j.s2*j.s2 + j.s3*j.s3); }

template<typename T, QType B>
T fabs (const Quaternion<T,B>& j)
{ return sqrt (norm(j)); }

// Return the positive definite square root of a Hermitian Quaternion
template<typename T, QType B>
Quaternion<T,B> sqrt (const Quaternion<T,B>& h)
{
  T root_det = sqrt( det(h) );
  T scalar = sqrt( 0.5 * (h.s0 + root_det) );

  if (scalar == 0.0)
    return Quaternion<T,B> (0.0);

  return Quaternion<T,B> (scalar, h.get_vector()/(2*scalar));
}

//! return a unitary matrix with rows equal to the eigenvectors of q
/*!
  Note that the eigenvalues of the hermitian input q=(I,Q,U,V)
  are I+p and I-p, where p is the magnitude of (Q,U,V)

  To derive the eigenvectors, start with

  For I+p: 
    e0=e1*(U-iV)/(p-Q)   [used when Q<0]
  or
    e1=e0*(U+iV)/(p+Q)   [used when Q>=0]

  For I-p:
    e0=-e1*(U-iV)/(p+Q)  [used when Q>=0]
  or
    e1=-e0*(U+iV)/(p-Q)  [used when Q<0]

  Plugging these directly into the columns of two Jones matrices (one
  for Q<0 and the other for Q>=0) and taking the Hermitian transpose
  produces the required basis transformation.  The determinant of the
  matrix in each case is

  det(J) = 2p/(p+Q)      [Q>=0]
  and
  det(J) = -2p/(p-Q)     [Q<0]

  Normalize each matrix by the sqrt of its determinant, and convert
  the result to a Unitary Quaternion ... done!
*/

template<typename T>
Quaternion<T,Unitary> eigen (const Quaternion<T,Hermitian>& q)
{
  T p = norm( q.get_vector() );

  /*
    q.s0 == 0 is a special case used by calculate_Jacobi and is required
    for the Jacobi method to work on Hermitian matrices.  Unfortunately,
    at the moment I can't remember how I derived calculate_Jacobi for
    Hermitian matrices and I can't update it to work in the case of 
    q.s1 < 0 ... another time ...

    Willem van Straten - 30 November 2009
  */
  if (q.s1 < 0 && q.s0 != 0)
  {
    T m = 1.0 / sqrt( 2.0*p*(p-q.s1) );
    return Quaternion<T,Unitary> (m*q.s3, -m*q.s2, -m*(p-q.s1), 0.0);
  }
  else
  {
    T m = 1.0 / sqrt( 2.0*p*(p+q.s1) );
    return Quaternion<T,Unitary> (m*(p+q.s1), 0.0, -m*q.s3, m*q.s2);
  }
}

//! Insertion operator
template<typename T>
std::ostream& operator<< (std::ostream& ostr, const Quaternion<T,Hermitian>& j)
{
  return ostr << "[h:" << j.s0 <<","<< j.s1 <<","<< j.s2 <<","<< j.s3 << "]";
}

//! Insertion operator
template<typename T>
std::ostream& operator<< (std::ostream& ostr, const Quaternion<T,Unitary>& j)
{
  return ostr << "[u:" << j.s0 <<","<< j.s1 <<","<< j.s2 <<","<< j.s3 << "]";
}

#endif  /* not __epsic_util_Quaternion_H defined */

