//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2011 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#ifndef __Dirac_H
#define __Dirac_H

#include "Pauli.h"
#include "Matrix.h"

namespace Dirac {

  typedef Matrix< 4,4,std::complex<double> > type;

  //! Get the specified basis matrix
  type matrix (unsigned i, unsigned j);

}


#endif
