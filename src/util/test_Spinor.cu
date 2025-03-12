//-*-C++-*-

/***************************************************************************
 *
 *   Copyright (C) 2025 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "Spinor.h"
#include <vector>

/*
  Multiply Spinors by Jones matrix
*/
__global__ void transform 
(
  const Jones<double>& jones,
  Spinor<double>* spinors,
  unsigned ndat
)
{
  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= ndat)
    return;

  spinors[idx] = jones * spinors[idx];
}

/*
  Compute array of instantaneous Stokes parameters from an array of spinors
*/
__global__ void compute_stokes 
(
  Vector<4,double>* stokes,
  Spinor<double>* spinors,
  unsigned ndat
)
{
  const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= ndat)
    return;

  compute_stokes(stokes[idx], spinors[idx]);
}

void launch ()
{
  Jones<double> jones;
  std::vector<Spinor<double>> spinors (1024);
  std::vector<Vector<4,double>> stokes (1024);

  dim3 threads (512);
  dim3 blocks (spinors.size() / threads.x);
  if (spinors.size() % threads.x)
    blocks.x ++;

  transform<<<blocks,threads>>> (jones, spinors.data(), spinors.size());
  compute_stokes<<<blocks,threads>>> (stokes.data(), spinors.data(), spinors.size());
}
