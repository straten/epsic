/***************************************************************************
 *
 *   Copyright (C) 2006 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "BoxMuller.h"
#include <math.h>
#include <stdlib.h>

BoxMuller::BoxMuller (long seed)
{
  have_one_ready = false;
  one_ready = 0;
  if (seed)
    srand48 (seed);
}

//! returns a random variable with a Gaussian distribution
float BoxMuller::evaluate ()
{
  if (have_one_ready) {
    /* use value from previous call */
    have_one_ready = false;
    return one_ready;
  }

  float v1, v2, w;

  do {
    v1 = 2.0 * drand48() - 1.0;
    v2 = 2.0 * drand48() - 1.0;
    w = v1 * v1 + v2 * v2;
  } while ( w >= 1.0 );
  
  w = sqrt( (-2.0 * log( w ) ) / w );

  have_one_ready = true;
  one_ready = v2 * w;

  return v1 * w;
  
}
