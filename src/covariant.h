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
private:
  covariant_mode* out[2];

  friend class covariant_mode;
  void get();

  // coefficient of mode intensity correlation
  double correlation;

protected:

  double get_correlation () const { return correlation; }

  //! Derived classes return a pair of mode intensities
  virtual void get_modulation (double& A, double& B) = 0;

public:

  //! Construct with correlation coefficient
  covariant_coordinator (double correlation);

  //! Virtual destructor (required for abstract base class)
  virtual ~covariant_coordinator () {}

  modulated_mode* get_modulated_mode (unsigned index, mode*);
};

class bivariate_lognormal_modes : public covariant_coordinator
{
  Matrix<2,2,double> meansq;
  Vector<2,double> mean;
  unsigned count;

  Matrix<2,2,double> correlator;
  void build ();
  bool built;
  double log_sigma[2];

  //! random number generator
  BoxMuller* normal;

protected:
  void get_modulation (double& A, double& B);

public:

  bivariate_lognormal_modes (double correlation) 
  : covariant_coordinator(correlation) 
  { built = false; set_beta(0,1.0); set_beta(1,1.0); }

  ~bivariate_lognormal_modes ();

  void set_beta (unsigned index, double);

  //! Return BoxMuller object used to generate normally distributed numbers
  virtual BoxMuller* get_normal () { return normal; }
  virtual void set_normal (BoxMuller* n) { normal = n; }

};

#endif
