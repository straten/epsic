//-*-C++-*-

/***************************************************************************
 *
 *   Copyright (C) 2025 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "Spinor.h"
#include "Jones.h"

#include <curand_kernel.h>

//! Functor to generate a random Spinor with known polarization state
/*! The random real and imaginary parts are normally distributed. */
class NormalSpinor
{
  unsigned long long base_seed = 0;
  unsigned long long sequence_index = 0;
  Jones<double> polarizer;

public:

  void set_seed (unsigned long long seed) { base_seed = seed; }
  void set_sequence_index (unsigned long long idx) { sequence_index = idx; }
  void set_polarizer (const Jones<double>& J) { polarizer = J; }

  __device__ Spinor<double> operator()();
};

__device__ Spinor<double> NormalSpinor::operator()()
{
  unsigned int local_tid = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned long long global_sequence_index = sequence_index + local_tid;

  // Initialize a cuRAND state for this specific element/thread
  curandState_t state;
  curand_init(base_seed, global_sequence_index, 0, &state);

  epsic::complex<double> x ( curand_normal_double(&state), curand_normal_double(&state) );
  epsic::complex<double> y ( curand_normal_double(&state), curand_normal_double(&state) );
  Spinor<double> e (x, y);
  return polarizer * e;
}