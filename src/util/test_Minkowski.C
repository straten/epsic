/***************************************************************************
 *
 *   Copyright (C) 2011 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "Minkowski.h"
#include "Dirac.h"
#include "Jacobi.h"

using namespace std;

int main () 
{
  // test that the Minkowski basis matrices are orthogonal to each other

  unsigned count = 0;

  for (unsigned i=0; i<4; i++)
   {
    for (unsigned j=i; j<4; j++)
    {
      Stokes<double> A = 0.0;
      Stokes<double> B = 0.0;

      A[i] = 1.0;
      B[j] = 1.0;

      Matrix<4,4,double> M1 = Minkowski::outer (A, B) + Minkowski::outer (B, A);
      cerr << count << "::" << M1 << endl;
      count ++;

      for (unsigned k=0; k<4; k++)
      {
        for (unsigned l=k; l<4; l++)
        {
          Stokes<double> A = 0.0;
          Stokes<double> B = 0.0;

          A[k] = 1.0;
          B[l] = 1.0;

          Matrix<4,4,double> M2 = Minkowski::outer (A, B) + Minkowski::outer (B, A);

          double tr = trace(M1 * M2);

          double expect = 0.0;
          if (i==k && j==l)
            expect = 2;
          if (i==j)
            expect *= 2;

          if (tr != expect)
          {   
            cerr << "fail " << i << " " << j << " " << k << " " << l << " " << tr << endl;
            return -1;
          }
        }
      }
    }
  }


  for (unsigned i=0; i<4; i++)
    for (unsigned j=0; j<4; j++)
    {
      Dirac::type M1 = Dirac::matrix (i, j);

      // test the mixed product property of the Kronecker product
      Dirac::type M2 = Dirac::matrix(i,0) * Dirac::matrix(0,j);
      if (M1 != M2)
      {
        cout << i << " " << j << endl << M1 << endl;
        cout << "!=\n" << M2 << endl << endl;
        return -1;
      }

      // test for symmetry
      bool symmetric = true;

      for (unsigned k=0; k<4 && symmetric; k++)
        for (unsigned l=k+1; l<4 && symmetric; l++)
          if ( M1[k][l] != M1[l][k] )
          {
            cerr << "not symmetric i=" << i << " j=" << j << endl;
            cerr << M1 << endl;
            symmetric = false;
          }
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
