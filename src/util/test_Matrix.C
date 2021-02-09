/***************************************************************************
 *
 *   Copyright (C) 2003 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
#include "Matrix.h"
#include "random.h"

using namespace std;

int main ()
{
  Matrix<3, 4, double> m1;
  Matrix<4, 5, double> m2;
  Vector<4, double> v;

  random_matrix (m1, 10.0);
  random_matrix (m2, 10.0);
  random_vector (v, 10.0);

  cerr << "m1=" << m1 << "\nm2=" << m2 << "\nm1*m2=" << m1*m2 << endl;

  cerr << "v=" << v << " m1*v=" << m1*v << endl;

  Matrix<3, 3, double> test;
  random_matrix (test, 10.0);

  Matrix<3, 3, double> inverse = inv (test);

  cerr << "M*inv(M)=" << test * inverse << endl;

  return 0;
}

