//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2006 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/BoxMuller.h

#ifndef __BoxMuller_h
#define __BoxMuller_h

//! Returns a random variable from a normal distribution
class BoxMuller {

 public:

  //! Default constructor
  BoxMuller (long seed = 0);

  //! returns a normal deviate with zero mean and unit variance
  float operator () () { return evaluate(); }

  float evaluate ();

 protected:

  //! one deviate is stored on each call to operator
  bool  have_one_ready;
  float one_ready;

};

#endif

