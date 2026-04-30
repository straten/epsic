//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2006 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/Convention.h

#ifndef __epsic_util_Convention_h
#define __epsic_util_Convention_h

#include <iostream>

//! Defines enumerations for various conventions relevant to polarimetry
/*! See Sections 5.1 Phase Convention and 5.4 Nominal Feed Configuration
    of [An Introduction to Single-antenna Radio Astronomical Polarimetry]
    (https://ui.adsabs.harvard.edu/abs/2026PASP..138a3001V/abstract)
    van Straten, W. (2026) PASP 138:013001 */
class Convention
{
public:

  //! The basis in which the electric field is represented
  /*! As described in 5.4.1. Feed Basis */
  enum Basis { Circular=0, Linear=1, Elliptical=2 };

  //! The hand of the basis
  /*! As described in 5.4.2. Feed Hand/Basis Reflection */
  enum Hand { Left=-1, Right=1 };

  //! The sign convention adopted for complex phase
  /*! As described in 5.1. Phase Convention */
  enum Argument { Conjugate=-1, Conventional=1 };
};

// //////////////////////////////////////////////////////////////////////

//! Basis output operator
std::ostream& operator << (std::ostream&, Convention::Basis);
//! Basis input operator
std::istream& operator >> (std::istream&, Convention::Basis&);

//! Hand output operator
std::ostream& operator << (std::ostream&, Convention::Hand);
//! Hand input operator
std::istream& operator >> (std::istream&, Convention::Hand&);

//! Argument output operator
std::ostream& operator << (std::ostream&, Convention::Argument);
//! Argument input operator
std::istream& operator >> (std::istream&, Convention::Argument&);

#endif

