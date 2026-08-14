/***************************************************************************
 *
 *   Copyright (C) 2006 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "BoxMuller.h"

constexpr float default_mean = 0.0f;
constexpr float default_stddev = 1.0f;

BoxMuller::BoxMuller (long seed)
: dist(default_mean, default_stddev)
{
  if (!seed)
  {
    std::random_device rd;
    seed = rd();
  }

  engine.seed(seed);
}

//! returns a random variable with a Gaussian distribution
float BoxMuller::evaluate ()
{
  return dist(engine);
}
