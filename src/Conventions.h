//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2006 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

/* $Source: /cvsroot/psrchive/psrchive/Util/units/Conventions.h,v $
   $Revision: 1.4 $
   $Date: 2006/10/06 21:13:55 $
   $Author: straten $ */

#ifndef __Conventions_h
#define __Conventions_h

#include <iostream>

namespace Signal {

  //! The basis in which the electric field is represented
  enum Basis { Circular=0, Linear=1, Elliptical=2 };

  //! The hand of the basis
  enum Hand { Left=-1, Right=1 };

  //! The complex phase of the basis
  enum Argument { Conjugate=-1, Conventional=1 };

  // //////////////////////////////////////////////////////////////////////

  //! Basis output operator
  std::ostream& operator << (std::ostream&, Basis);
  //! Basis input operator
  std::istream& operator >> (std::istream&, Basis&);

  //! Hand output operator
  std::ostream& operator << (std::ostream&, Hand);
  //! Hand input operator
  std::istream& operator >> (std::istream&, Hand&);

  //! Argument output operator
  std::ostream& operator << (std::ostream&, Argument);
  //! Argument input operator
  std::istream& operator >> (std::istream&, Argument&);

}

#endif

