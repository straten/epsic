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

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/transform.h>

#include <iostream>
#include <ctime>

using namespace std;

int test_PolarizedNoise(PolarizedNoise<double>& generator)
{
  size_t nsamp = 1024*1024;

  {
    cerr << "generating on device using thrust" << endl;

    PolarizedNoise<float> fgenerator = generator;

    thrust::device_vector<Spinor<float>> e_field (nsamp);
    thrust::generate(e_field.begin(), e_field.end(), fgenerator);

    cerr << "transforming on device using thrust" << endl;

    thrust::device_vector<Stokes<float>> stokes (nsamp);
    thrust::transform(e_field.begin(), e_field.end(), stokes.begin(), instantaneous_stokes<float>());

    cerr << "summing on device using thrust" << endl;

    Stokes<float> sum = thrust::reduce(stokes.begin(), stokes.end());

    cerr << "thrust device mean=" << sum/nsamp << endl << endl;
  }

  cerr << "generating on host using thrust" << endl;

  thrust::host_vector<Spinor<double>> e_field (nsamp);
  thrust::generate(e_field.begin(), e_field.end(), generator);

  cerr << "transforming on host using thrust" << endl;

  thrust::host_vector<Stokes<double>> stokes (nsamp);
  thrust::transform(e_field.begin(), e_field.end(), stokes.begin(), instantaneous_stokes<double>());

  cerr << "summing on host using thrust" << endl;

  Stokes<double> sum = thrust::reduce(stokes.begin(), stokes.end());

  cerr << "thrust host mean=" << sum/nsamp << endl;

  return 0;
}
