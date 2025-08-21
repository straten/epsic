//-*-C++-*-

/***************************************************************************
 *
 *   Copyright (C) 2025 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "PolarizedNoise.h"
#include "Stokes.h"
#include "Quaternion.h"
#include "Pauli.h"

#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>

using namespace std;

// defined in test_PolarizedNoise.cu
extern int test_PolarizedNoise(PolarizedNoise<double>& generator);

int main()
{
  PolarizedNoise<double> generator;
  generator.set_seed(time(NULL));

  Stokes<double> mean (1.0, -0.1, 0.2, 0.6);

  cerr << "Expected Stokes parameters: " << mean << endl << endl;

  Quaternion<double,Hermitian> root = 0.5 * sqrt (natural(mean));
  generator.set_polarizer(convert(root));

  size_t nsamp = 1024*1024;

  cerr << "generating on host using STL" << endl;

  std::vector<Spinor<double>> e_field (nsamp);
  std::generate(e_field.begin(), e_field.end(), generator);

  cerr << "transforming on host using STL" << endl;

  std::vector<Stokes<double>> stokes (nsamp);
  std::transform(e_field.begin(), e_field.end(), stokes.begin(), instantaneous_stokes<double>());

  cerr << "summing on host using STL" << endl;

  Stokes<double> sum = std::reduce(stokes.begin(), stokes.end());

  cerr << "STL host mean=" << sum/nsamp << endl << endl;

  // two more times using thrust
  return test_PolarizedNoise(generator);
}
