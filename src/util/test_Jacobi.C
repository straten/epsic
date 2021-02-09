/***************************************************************************
 *
 *   Copyright (C) 2004 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "Jacobi.h"
#include "MatrixTest.h"

double real (double x) { return x; }
double conj (double x) { return x; }

// specialize for complex matrices
template <unsigned RC, typename T>
void test_Jacobi (Matrix<RC,RC,T>& matrix, double tolerance)
{
  unsigned i,j;

  // make the incoming matrix Hermitian
  for (i=0; i<RC; i++) {
    matrix[i][i] = real(matrix[i][i]);
    for (j=i+1; j<RC; j++)
      matrix[i][j] = conj(matrix[j][i]);
  }

  // make a copy
  Matrix<RC, RC, T> temp = matrix;

  Matrix<RC, RC, T>  eigenvectors;
  Vector<RC, double> eigenvalues;

  Jacobi (temp, eigenvectors, eigenvalues);
  temp = eigenvectors * matrix * herm(eigenvectors);

  // ensure that the result is diagonal and that the diagonal elements
  // are equal to the eigenvalues

  for (i=0; i<RC; i++)
    for (j=0; j<RC; j++) {

      double should_be = 0;
      if (i == j)
	should_be = eigenvalues[i];

      if (norm(temp[i][j] - should_be) > tolerance) {
	std::cerr << "test_Jacobi t[" << i << "][" << j << "]=" << temp[i][j]
		  << " != " << should_be << std::endl;

	std::cerr << "P = " << eigenvectors << std::endl;
	std::cerr << "P^T A P = " << temp << std::endl;

	throw std::string ("test_Jacobi error");
      }
    }

}

template<typename T, unsigned dim> 
void runtest (unsigned loops, double tolerance)
{
  std::cerr << loops << " " << dim << "x" << dim << " Jacobi solutions" << std::endl;

  for (unsigned iloop=0; iloop<loops; iloop++) {

    Matrix<dim,dim,T> m1;
    random_matrix (m1, 10.0);

    test_Jacobi (m1, tolerance);

  }
}

int main () try {

  unsigned loops = 100000;
  double tolerance = 1e-12;

  std::cerr << "Testing symmetric: ";
  runtest<double,7> (loops, tolerance);

  std::cerr << "Testing Hermitian: ";
  runtest<std::complex<double>,5> (loops, tolerance);

  std::cerr << "Successful tests" << std::endl;
  return 0;

}
catch (std::string& error) {
  std::cerr << "test_Jacobi error "<< error << std::endl;
  return -1;
}
