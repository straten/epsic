#include "Pauli.h"
#include "random.h"

#include <iostream>
using namespace std;

int main ()
{
  Vector< 4, complex<double> > input1;
  Vector< 4, complex<double> > input2;

  for (unsigned i=0; i < 1000; i++) {

    // simulate some complex Stokes parameters to enter the inner product
    random_vector (input1, 10.0);
    random_vector (input2, 10.0);

    // Stokes 4-vectors
    Stokes< complex<double> > S1 = input1;
    Stokes< complex<double> > S2 = input2;

    // equivalent coherency matrices
    Jones<double> J1 = convert(S1);
    Jones<double> J2 = convert(S2);

    double inner1 = (S1 * conj(S2)).real();
    double inner2 = 2 * trace(J1*herm(J2)).real();

    double diff = fabs( (inner1 - inner2)/inner1 );

    if (diff > 1e-12) {
      cerr << "inner1=" << inner1 << endl
	   << "inner2=" << inner2 << endl
	   << "diff=" << diff << endl
	   << "S1=" << S1 << endl
	   << "S2=" << S2 << endl
	   << "conj(S2)=" << conj(S2) << endl
	   << "J1=" << J1 << endl
	   << "J2=" << J2 << endl
	   << "herm(J2)=" << herm(J2) << endl;
      return -1;
    }

  }

  cerr << "inner product pass: Stokes parameters and Jones matrices" << endl;

  return 0;
}

