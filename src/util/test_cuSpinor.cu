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

#include <iostream>
#include <ctime>

int test_cuSpinor()
{
  size_t nsamp = 1024*1024;
  thrust::device_vector<Spinor<double>> e_field (nsamp);

  NormalSpinor generator;
  generator.set_seed(time(NULL));

  Stokes<double> mean (1.0, 0.0, 0.2, 0.7);
  Quaternion<double,Hermitian> root = sqrt (natural(mean));
  generator.set_polarizer(convert(root));

  thrust::generate(e_field.begin(), e_field.end(), generator);

  thrust::host_vector<Spinor<double>> on_host = e_field;

  for (int i = 0; i < 10; ++i)
  {
    std::cout << "  " << i << ": " << on_host[i] << std::endl;
  }
  std::cout << std::endl;

  return 0;
}
