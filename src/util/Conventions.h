//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2006 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/Conventions.h

#ifndef __epsic_util_Conventions_h
#define __epsic_util_Conventions_h

#include <iostream>

class Signal
{
public:

  //! The basis in which the electric field is represented
  enum Basis { Circular=0, Linear=1, Elliptical=2 };

  //! The hand of the basis
  enum Hand { Left=-1, Right=1 };

  //! The complex phase of the basis
  enum Argument { Conjugate=-1, Conventional=1 };
};

// //////////////////////////////////////////////////////////////////////

//! Basis output operator
std::ostream& operator << (std::ostream&, Signal::Basis);
//! Basis input operator
std::istream& operator >> (std::istream&, Signal::Basis&);

//! Hand output operator
std::ostream& operator << (std::ostream&, Signal::Hand);
//! Hand input operator
std::istream& operator >> (std::istream&, Signal::Hand&);

//! Argument output operator
std::ostream& operator << (std::ostream&, Signal::Argument);
//! Argument input operator
std::istream& operator >> (std::istream&, Signal::Argument&);

#endif

