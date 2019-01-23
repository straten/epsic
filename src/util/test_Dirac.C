/***************************************************************************
 *
 *   Copyright (C) 2011 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "Dirac.h"
#include "Jacobi.h"

using namespace std;

int main () 
{
  for (unsigned i=0; i<4; i++)
    for (unsigned j=0; j<4; j++)
    {
      Dirac::type M = Dirac::matrix (i, j);
      cout << i << " " << j << endl << M << endl;
      cout << "==\n" << Dirac::matrix(i,0)*Dirac::matrix(0,j) << endl << endl;
    }

  for (unsigned i=0; i<2; i++)
  {
    for (unsigned j=1; j<4; j++)
    {
      Dirac::type covar (0);
      covar[0][0] = covar[j][j] = 1.0;

      if (i == 0)
	covar[j][0] = covar[0][j] = 1.0;

    cout << "covar=\n" << covar << endl;

    Dirac::type N (0);
    for (unsigned i=0; i<4; i++)
      for (unsigned j=0; j<4; j++)
	N += Dirac::type(covar[i][j]) * Dirac::matrix (i, j);

    cout << "N=\n" << N << endl << endl;

    Dirac::type evec;
    Vector<4,double> eval;
    Jacobi (N, evec, eval);

    for (unsigned i=0; i<4; i++)
      cout << eval[i] << " " << evec[i] << endl;
    }
  }

  return 0;
}
