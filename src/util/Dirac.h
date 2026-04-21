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

//! Generates Dirac matrices as Kronecker products of Pauli+identity matrices
class Dirac 
{
public:
  typedef Matrix< 4,4,std::complex<double> > type;

  //! Get the specified Dirac basis matrix
  /*! The result is the Kronecker product: $\sigma_i \ocross \sigma_j$*/
  static type matrix (unsigned i, unsigned j);
};

#endif
