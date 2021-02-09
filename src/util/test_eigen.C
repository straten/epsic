/***************************************************************************
 *
 *   Copyright (C) 2004 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "MatrixTest.h"
#include "Pauli.h"

using namespace std;

template <typename T> T sqr (T x) { return x*x; }

template <typename T>
void test_eigen(const Quaternion<T, Hermitian>& q, double tolerance)
{
  // cerr << "test_eigen: q=" << q << endl;

  Quaternion<T, Unitary> eq = eigen (q);

  // cerr << "test_eigen: eq=" << eq << endl;

  double determinant = det(eq);
  if ( sqr(determinant - 1.0) > tolerance ) {
    cerr << "test_eigen det(unitary)=" << determinant <<" != 1"<<endl;
    throw string ("test_eigen error");
  }

  Jones<T> R   = convert(eq);
  Jones<T> rho = convert(q);

  // cerr << "test_eigen: rho=" << rho << endl;

  Jones<T> result = R * rho * herm(R);

  // cerr << "test_eigen: res=" << result << endl;

  double off = norm(result(0,1)) + norm(result(1,0));

  if ( sqrt(off) > tolerance )
  {
    cerr << "tol=" << tolerance << endl;
    cerr << "off=" << sqrt(off) << endl;

    cerr << "test_eigen result not diagonal" <<endl;
    throw string ("test_eigen error");
  }

  double det1 = det(q);
  double det2 = result(0,0).real() * result(1,1).real();

  if ( sqr(det1 - det2) > tolerance )
  {
    cerr << "tol=" << tolerance << endl;
    cerr << "off=" << sqr(det1 - det2) << endl;

    cerr << "test_eigen det(q)=" << det1 <<" != l1*l2=" << det2 <<endl;
    throw string ("test_eigen error");
  }

}

// specialize the test_matrix template function defined in MatrixTest.h
template <typename T>
void test_matrix (const Quaternion<T, Hermitian>& q1, 
		  const Quaternion<T, Hermitian>& q2,
		  const double& c0, const double& c1,
		  double scale, double tolerance, bool verbose)
{
  if (verbose)
    cerr << "test_matrix eigen quaternion" << endl;

  tolerance *= 500;

  test_eigen (q1, tolerance);
  test_eigen (q2, tolerance);

}


int main () 
{
  unsigned loops = 1024 * 1024;

  MatrixTest <Quaternion<double,Hermitian>, double, double > test;

  try {
    cerr << "Testing " << loops << " eigen solutions" << endl;
    test.runtest (loops);
  }
  catch (string& error) {
    cerr << error << endl;
    return -1;
  }

  return 0;
}

