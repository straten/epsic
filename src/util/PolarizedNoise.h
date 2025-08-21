//-*-C++-*-

/***************************************************************************
 *
 *   Copyright (C) 2025 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/PolarizedNoise.h

#ifndef _epsic_PolarizedNoise_H
#define _epsic_PolarizedNoise_H

#include "Spinor.h"
#include "Jones.h"

#include <random>
#include <iostream>

#ifdef __CUDACC__
#include <curand_kernel.h>
#else
#define __host__
#endif

//! Functor generates a random Spinor with known polarization state
/*! The random real and imaginary parts are normally distributed. */
template<typename T>
class PolarizedNoise
{
  unsigned long long base_seed = 0;
  unsigned long long sequence_index = 0;
  Jones<T> polarizer;

  //! Random number generator used on the host
  mutable std::mt19937* host_generator = nullptr;

public:

  //! Host constructor intializes the host generator
  __host__ PolarizedNoise ()
  {
    // Initialize the host-side generator
    std::random_device rd;
    host_generator = new std::mt19937(rd());
  }

  //! Generate a polarized random Spinor
  __prefix__ Spinor<T> operator()();

  //! Set the Jones matrix used to polarize the generated Spinor instances
  void set_polarizer (const Jones<T>& J) { polarizer = J; }

  void set_seed (unsigned long long seed) 
  {
    base_seed = seed;
    if (host_generator)
      host_generator->seed(seed);
  }

  void set_sequence_index (unsigned long long idx) { sequence_index = idx; }
};


template<typename T>
__prefix__ Spinor<T> PolarizedNoise<T>::operator()()
{
#ifdef __CUDA_ARCH__
  unsigned int local_tid = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned long long global_sequence_index = sequence_index + local_tid;

  curandState_t state;
  curand_init(base_seed, global_sequence_index, 0, &state);

  epsic::complex<double> x ( curand_normal_double(&state), curand_normal_double(&state) );
  epsic::complex<double> y ( curand_normal_double(&state), curand_normal_double(&state) );
#else
  std::normal_distribution<T> normal;
  epsic::complex<T> x ( normal(*host_generator), normal(*host_generator) );
  epsic::complex<T> y ( normal(*host_generator), normal(*host_generator) );
#endif

  Spinor<T> e (x, y);
  return polarizer * e;
}


#endif // _epsic_PolarizedNoise_H