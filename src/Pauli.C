/***************************************************************************
 *
 *   Copyright (C) 2004 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
#include "Pauli.h"

Basis<double>& Pauli::basis()
{
  /*
    The order of construction of global variables at runtime is not
    guaranteed, but the following pointer is initialized at compile time
  */

  static Basis<double>* basis = 0;
  if (!basis)
    basis = new Basis<double>;

  return *basis;
}

Jones<double> Pauli::matrix (unsigned i)
{
  Quaternion<double,Hermitian> q;
  q[i] = 1;
  return convert (q);
}
