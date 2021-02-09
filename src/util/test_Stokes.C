/***************************************************************************
 *
 *   Copyright (C) 2003 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
#include "Pauli.h"

#include <iostream>
using namespace std;

int main ()
{
  Stokes< complex< Estimate<double> > > test_default_constructor;

  Jones<double> jones;
  Stokes<double> in;

  Stokes<double> out = transform (in, jones);

  random_value (in, 10.0, 0.99);
  if (in.abs_vect() > 9.90) {
    cerr << "error: random_value" << endl;
    return -1;
  }

  Stokes<double> S_test (5, 1, 2, 3);
  Jones<double> J_test = convert(S_test);

  complex<double> half (0.5);

  if (J_test(0,0) != complex<double> (5+1)*half) {
    cerr << "error: I+Q" << endl;
    return -1;
  }
  if (J_test(1,1) != complex<double> (5-1)*half) {
    cerr << "error: I-Q" << endl;
    return -1;
  }
  if (J_test(1,0) != complex<double> (2,3)*half) {
    cerr << "error: U,V" << endl;
    return -1;
  }
  if (J_test(0,1) != complex<double> (2,-3)*half) {
    cerr << "error: U,-V" << endl;
    return -1;
  }

  Stokes<double> S_back = coherency (J_test);
  if (S_back != S_test) {
    cerr << "error: back=" << S_back 
	 << "   !=  test=" << S_test << endl;
    return -1;
  }

  Pauli::basis().set_basis( Signal::Circular );

  J_test = convert(S_test);

  if (J_test(0,0) != complex<double> (5+3)*half) {
    cerr << "error: I+V" << endl;
    return -1;
  }
  if (J_test(1,1) != complex<double> (5-3)*half) {
    cerr << "error: I-V" << endl;
    return -1;
  }
  if (J_test(1,0) != complex<double> (1,2)*half) {
    cerr << "error: Q,U" << endl;
    return -1;
  }
  if (J_test(0,1) != complex<double> (1,-2)*half) {
    cerr << "error: Q,-U" << endl;
    return -1;
  }

  cerr << "Stokes class passes tests" << endl;

  return 0;
}

