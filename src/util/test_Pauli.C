/***************************************************************************
 *
 *   Copyright (C) 2003 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
#include "MatrixTest.h"
#include "Pauli.h"

using namespace std;

template <typename T> T sqr (T x) { return x*x; }


template <typename T, QBasis B>
void test_det (const Quaternion<T, B>& q, double tolerance)
{
  Jones<T> j = convert(q);

  T det1 = det(j).real();
  T det2 = det(q);

  if ( sqr(det1 - det2) > tolerance ) {
    cerr << "test_det det(Jones)=" << det1 << " != det(Quaternion)=" << det2
	 << endl;
    throw string ("test_Pauli determinant error");
  }
}

template <typename T>
void test_polar (const Jones<T>& j, double tolerance)
{
  // pull the jones matrix apart into its Hermitian and Unitary components
  std::complex<T> d;
  Quaternion<T, Hermitian> hq;
  Quaternion<T, Unitary> uq;

  polar (d, hq, uq, j);

  double determinant = 0;

  determinant = det(hq);
  if ( sqr(determinant - 1.0) > tolerance ) {
    cerr << "test_matrix polar det(hermitian)=" << determinant <<" != 1"<<endl;
    throw string ("test_matrix polar decomposition error");
  }

  determinant = det(uq);
  if ( sqr(determinant - 1.0) > tolerance ) {
    cerr << "test_matrix polar det(unitary)=" << determinant <<" != 1"<<endl;
    throw string ("test_matrix polar decomposition error");
  }


  // reconstruct the Jones matrix
  Jones<T> oj = d * (hq * uq);

  // test that you get the same matrix back
  double diff = norm(oj-j)/norm(oj);
  if ( diff > tolerance ) {
    cerr << "test_matrix polar\n"
      " out=" << oj << " !=\n"
      "  in=" << j << endl << " diff=" << diff << endl;
    throw string ("test_matrix polar decomposition error");
  }

  double value;
  random_value (value, 10.0);
  
  test_det (hq*value, tolerance);
  test_det (uq*value, tolerance);
}

// specialize the test_matrix template function defined in MatrixTest.h
template <typename T>
void test_matrix (const Jones<T>& j1, const Jones<T>& j2,
		  const Quaternion<T>& q, const std::complex<T>& c,
		  double scale, double tolerance, bool verbose)
{
  if (verbose)
    cerr << "test_matrix Jones/Quaternion Isomorphism" << endl;

  test_polar (j1, tolerance);
  test_polar (j2, tolerance);
}


int main () 
{
  complex<double> i (0, 1);

  Jones<double> sigma_1 (1, 0,
                         0, -1);

  Jones<double> sigma_2 (0, 1,
                         1, 0);

  Jones<double> sigma_3 (0, -i,
                         i, 0);

  //
  // test Unitary transformations
  //
  Quaternion<double,Unitary> R_1 (0, 1, 0, 0);
  Jones<double> JR_1 = convert (R_1);
  cerr << "R_1=" << JR_1 << endl;
  if (JR_1 != i * sigma_1)
  {
    cerr << "R_1 does not equal expected value\n" << i*sigma_1 << endl;
    return -1;
  }

  Quaternion<double,Unitary> R_2 (0, 0, 1, 0);
  Jones<double> JR_2 = convert (R_2);
  cerr << "R_2=" << JR_2 << endl;
  if (JR_2 != i * sigma_2)
  {
    cerr << "R_2 does not equal expected value\n" << i*sigma_2 << endl;
    return -1;
  }

  Quaternion<double,Unitary> R_3 (0, 0, 0, 1);
  Jones<double> JR_3 = convert (R_3);
  cerr << "R_3=" << JR_3 << endl;
  if (JR_3 != i * sigma_3)
  {
    cerr << "R_3 does not equal expected value\n" << i*sigma_3 << endl;
    return -1;
  }

  //
  // test Hermitian transformations
  //
  Jones<double> JB_1 = Pauli::matrix(1);
  cerr << "B_1=" << JB_1 << endl;
  if (JB_1 != sigma_1)
  {
    cerr << "B_1 does not equal expected value\n" << sigma_1 << endl;
    return -1;
  }

  Jones<double> JB_2 = Pauli::matrix(2);
  cerr << "B_2=" << JB_2 << endl;
  if (JB_2 != sigma_2)
  {
    cerr << "B_2 does not equal expected value\n" << sigma_2 << endl;
    return -1;
  }

  Jones<double> JB_3 = Pauli::matrix(3);
  cerr << "B_3=" << JB_3 << endl;
  if (JB_3 != sigma_3)
  {
    cerr << "B_3 does not equal expected value\n" << sigma_3 << endl;
    return -1;
  }


  unsigned loops = 1024 * 1024;

  MatrixTest <Jones<double>, Quaternion<double>, std::complex<double> > test;

  try {
    cerr << "Testing " << loops << " Pauli conversion variations" << endl;
    test.runtest (loops);
  }
  catch (string& error) {
    cerr << error << endl;
    return -1;
  }

  return 0;
}
