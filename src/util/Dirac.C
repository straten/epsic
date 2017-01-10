/***************************************************************************
 *
 *   Copyright (C) 2011 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "Dirac.h"

Dirac::type Dirac::matrix (unsigned i, unsigned j)
{
  typedef Matrix< 2,2,std::complex<double> > pauli_type;
  pauli_type sigma_i = Pauli::matrix(i);
  pauli_type sigma_j = Pauli::matrix(j);

  return direct (sigma_i, sigma_j);
}
