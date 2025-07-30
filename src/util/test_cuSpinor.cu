//-*-C++-*-

/***************************************************************************
 *
 *   Copyright (C) 2025 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "cuSpinor.h"
#include "Stokes.h"
#include "Quaternion.h"
#include "Pauli.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/transform.h>

#include <iostream>
#include <ctime>

template<typename T>
class instantaneous_stokes
{
  __prefix__ Stokes<T> operator() (const Spinor<T>& e)
  {
    Vector<4,T> tmp;
    compute_stokes (tmp, e);
    return tmp;
  }
};

int test_cuSpinor()
{
  NormalSpinor generator;
  generator.set_seed(time(NULL));

  Stokes<double> mean (1.0, 0.0, 0.2, 0.7);
  Quaternion<double,Hermitian> root = 0.5 * sqrt (natural(mean));
  generator.set_polarizer(convert(root));

  size_t nsamp = 1024*1024;
  thrust::device_vector<Spinor<double>> e_field (nsamp);
  thrust::generate(e_field.begin(), e_field.end(), generator);

  thrust::device_vector<Stokes<double>> stokes (nsamp);
  thrust::transform(e_field.begin(), e_field.end(), stokes.begin(), instantaneous_stokes<double>());

  Stokes<double> sum = thrust::reduce(stokes.begin(), stokes.end());

  std::cout << "mean=" << sum/nsamp << std::endl;

  return 0;
}
