/***************************************************************************
 *
 *   Copyright (C) 2003 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "Quaternion.h"

using namespace std;

#include "MatrixTest.h"

int main () 
{
  unsigned loops = 1024 * 1024;

  MatrixTest <Quaternion<double,Unitary>,
    Quaternion<double,Unitary>, double> testu;

  try {
    cerr << "Testing " << loops << " Unitary Quaternion variations" << endl;
    testu.runtest (loops);
  }
  catch (string& error) {
    cerr << error << endl;
    return -1;
  }




  MatrixTest <Quaternion<std::complex<double>,Hermitian>, 
    Quaternion<std::complex<double>,Hermitian>, std::complex<double> > testh;

  try {
    cerr 
      << "Testing " << loops << " Hermitian Biquaternion variations" << endl;
    testh.runtest (loops);
  }
  catch (string& error) {
    cerr << error << endl;
    return -1;
  }


  MatrixTest <Quaternion<std::complex<double>, Unitary>, 
    Quaternion<std::complex<double>, Unitary>, std::complex<double> > testub;

  try {
    cerr 
      << "Testing " << loops << " Unitary Biquaternion variations" << endl;
    testub.runtest (loops);
  }
  catch (string& error) {
    cerr << error << endl;
    return -1;
  }

  return 0;
}
