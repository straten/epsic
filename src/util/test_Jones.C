/***************************************************************************
 *
 *   Copyright (C) 2003 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
#include "MatrixTest.h"
#include "Jones.h"

using namespace std;

int main () 
{
  unsigned loops = 1024 * 1024;

  MatrixTest <Jones<double>, Jones<double>, std::complex<double> > test;

  try {
    cerr << "Testing " << loops << " Jones matrix variations" << endl;
    test.runtest (loops);
  }
  catch (string& error) {
    cerr << error << endl;
    return -1;
  }

  return 0;
}
