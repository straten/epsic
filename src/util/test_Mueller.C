#include "Pauli.h"

#include <iostream>

using namespace std;

int main ()
{
  for (unsigned i=0; i < 1000; i++) {

    Jones<double> J;
    random_vector (J, 10.0);

    Stokes<double> input;

    random_vector (input, 10.0);

    Stokes<double> output1 = transform (input, J);

    Matrix<4,4,double> M = Mueller (J);

    Stokes<double> output2 = M * input;

    Stokes<double> diff = output2 - output1;

    if (norm(diff)/norm(output1) > 1e-14) {
      cerr << "Error: out1=" << output1 << endl
	   << "    != out2=" << output2 << endl;
      return -1;
    }
  }

  cerr << "Mueller passes test" << endl;

  return 0;
}
