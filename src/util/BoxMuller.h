//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2006-2026 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/BoxMuller.h

#ifndef __epsic_util_BoxMuller_h
#define __epsic_util_BoxMuller_h

#include <random>

//! Returns a random variable from a normal distribution
/*! Uses std::normal_distribution, which typically implements the Marsaglia polar method,
    an optimized variant of Box-Muller that uses rejection sampling within a unit circle
    to avoid sine and cosine calculations, requiring only std::log and std::sqrt.
    The previous implementation of this class also used the Marsaglia polar method;
    however, it also used drand48 to generate 64-bit double-precision floating point
    values, and only 32-bit floats are returned.
  */
class BoxMuller
{
  // Mersenne Twister random number engine
  std::mt19937 engine;

  // Uniform distribution, constrained to output floats
  std::normal_distribution<float> dist;

  public:

  //! Default constructor
  BoxMuller (long seed = 0);

  //! returns a normal deviate with zero mean and unit variance
  float operator () () { return evaluate(); }

  //! returns a normal deviate with zero mean and unit variance
  float evaluate ();
};

#endif

