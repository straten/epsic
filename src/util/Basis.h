//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2004 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/Basis.h

#ifndef __Basis_H
#define __Basis_H

#include "Matrix.h"
#include "Conventions.h"

//! Defines the basis in which the electric field is represented
template<typename T>
class Basis {

public:

  //! Default constructor
  Basis () { set_basis (Signal::Linear); }

  //! set basis to circular or linear
  void set_basis (Signal::Basis basis);

  //! get the current basis
  Signal::Basis get_basis () const { return basis; }

  //! set basis to elliptical
  void set_basis (double orientation, double ellipticity);

  //! get the orientation
  double get_orientation () const { return orientation; }

  //! get the ellipticity
  double get_ellipticity () const { return ellipticity; }

  //! get the basis vector
  Vector<3,T> get_basis_vector (unsigned iaxis) const { return into[iaxis]; }

  //! convert the input vector into the basis
  template<typename U>
  Vector<3,U> get_in (const Vector<3,U>& vect) const { return into * vect; }

  //! convert the input vector out of the basis
  template<typename U>
  Vector<3,U> get_out (const Vector<3,U>& vect) const { return outof * vect; }

protected:

  //! The basis code
  Signal::Basis basis;

  //! The orientation 
  double orientation;

  //! The ellipticity
  double ellipticity;

  //! The basis matrix
  Matrix<3, 3, T> into;

  //! The transpose of the basis matrix
  Matrix<3, 3, T> outof;
};

/*! Given the orientation and ellipticity in radians, calculate the basis
  vectors used to convert Stokes parameters to coherency matrix. */
template<typename T>
void Basis<T>::set_basis (double _orientation, double _ellipticity)
{
  orientation = _orientation;
  ellipticity = _ellipticity;

  double cos_2o = cos (2.0*orientation);
  double sin_2o = sin (2.0*orientation);
  double cos_2e = cos (2.0*ellipticity);
  double sin_2e = sin (2.0*ellipticity);

  into[0] = Vector<3,T> (cos_2o*cos_2e, sin_2o, -cos_2o*sin_2e);
  into[1] = Vector<3,T> (-sin_2o*cos_2e, cos_2o, sin_2o*sin_2e);
  into[2] = Vector<3,T> (sin_2e, 0, cos_2e);

  basis = Signal::Elliptical;

  outof = transpose (into);
}

//! set basis to circular or linear
template<typename T>
void Basis<T>::set_basis (Signal::Basis _basis)
{
  basis = _basis;

  switch (basis)  {

  case Signal::Linear:
    into[0] = Vector<3,T>::basis (0); // hat q
    into[1] = Vector<3,T>::basis (1); // hat u
    into[2] = Vector<3,T>::basis (2); // hat v
    orientation = 0;
    ellipticity = 0;
    break;

  case Signal::Circular:
    into[0] = Vector<3,T>::basis (1); // hat q
    into[1] = Vector<3,T>::basis (2); // hat u
    into[2] = Vector<3,T>::basis (0); // hat v
    orientation = 0.25*M_PI;
    ellipticity = 0.25*M_PI;
    break;

  default:
    throw std::runtime_error ("Basis::set_basis unrecognized basis");

  }

  outof = transpose (into);
}

#endif

