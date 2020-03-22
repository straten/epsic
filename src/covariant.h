//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2020 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/covariant.h

#ifndef __covariant_H
#define __covariant_H

#include "modulated.h"
#include "Matrix.h"

#include <queue>

/***************************************************************************
 *
 *  models covariant mode intensities
 *
 ***************************************************************************/

class covariant_coordinator;

class covariant_mode : public modulated_mode
{
  friend class covariant_coordinator;

  covariant_coordinator* coordinator;
  std::queue<double> amps;

  double mean;
  double variance;

public:

  // initialize based on the modulation index beta
  covariant_mode (mode* s) : modulated_mode (s) { coordinator = 0; }

  // return a random scalar modulation factor
  double modulation ();

  double get_mod_mean () const { return mean; }
  
  double get_mod_variance () const { return variance; }

};

class covariant_coordinator
{
  modulated_mode* a_in;
  covariant_mode* a_out;

  modulated_mode* b_in;
  covariant_mode* b_out;

  friend class covariant_mode;
  void get();

  Matrix<2,2,double> covar;

public:

  covariant_coordinator (double covariance);

  void set_modeA_input (modulated_mode*);
  void set_modeB_input (modulated_mode*);

  modulated_mode* get_modeA_output ();
  modulated_mode* get_modeB_output ();
};

#endif
