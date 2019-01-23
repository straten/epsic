/***************************************************************************
 *
 *   Copyright (C) 2004 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
#include "Basis.h"
#include <iostream>

using namespace std;

int main ()
{
  Basis<double> basis1;
  basis1.set_basis (0.25*M_PI, 0.25*M_PI);

  Basis<double> basis2;
  basis2.set_basis (Signal::Circular);

  for (unsigned i=0; i<3; i++)  {
    double dot = basis1.get_basis_vector(i) * basis2.get_basis_vector(i);
    if (dot != 1) {
      cerr << "basis1 != basis2" << endl;
      return -1;
    }
  }

  cerr << "Basis::Circular basis has expected orientation=ellipticity=pi/4"
       << endl;

  return 0;
}
