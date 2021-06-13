//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2016 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/sample.h

#ifndef __sample_H
#define __sample_H

#include "mode.h"

#include <cstdlib>
#include <vector>

/***************************************************************************
 *
 *  a Stokes sample
 *
 ***************************************************************************/

class sample
{
public:

  unsigned sample_size;

  sample() { sample_size = 1; }

  virtual ~sample () {}

  virtual Stokes<double> get_Stokes () = 0;
  virtual Vector<4, double> get_mean () = 0;
  virtual Matrix<4,4, double> get_covariance () = 0;
  virtual Matrix<4,4, double> get_crosscovariance (unsigned ilag) = 0;

  Matrix<4,4, double> get_covariance (mode* s, unsigned sample_size);

  Matrix<4,4, double> get_crosscovariance (mode* s, unsigned at_lag,
					   unsigned sample_size);
};

/***************************************************************************
 *
 *  a single source of electromagnetic radiation
 *
 ***************************************************************************/

class single : public sample
{
public:
  mode* source;

  single (mode* s) { source = s; }

  ~single () { delete source; }

  virtual Stokes<double> get_Stokes_instance ()
  {
    Spinor<double> e = source->get_field();
    Vector<4, double> tmp;
    compute_stokes (tmp, e);
    return tmp;
  }

  Stokes<double> get_Stokes ()
  {
    Stokes<double> result;
    for (unsigned i=0; i<sample_size; i++)
      result += get_Stokes_instance();
    result /= sample_size;
    return result;
  }

  Vector<4, double> get_mean ()
  {
    return source->get_mean();
  }

  Matrix<4,4, double> get_covariance ()
  {
    return sample::get_covariance (source, sample_size);
  }

  Matrix<4,4, double> get_crosscovariance (unsigned ilag)
  {
    return sample::get_crosscovariance (source, ilag, sample_size);
  }

};


/***************************************************************************
 *
 *  a combination of two sources of electromagnetic radiation
 *
 ***************************************************************************/

class combination : public sample
{
protected:
  double intensity_covariance;

public:
  mode* A;
  mode* B;

  combination () { A = new mode; B = new mode; intensity_covariance = 0; }

  virtual void set_normal (BoxMuller* n)
  { A->set_normal(n); B->set_normal(n); }

  //! Covariant mode intensities
  void set_intensity_covariance (double covar) { intensity_covariance = covar; }

  Matrix<4,4, double> get_crosscovariance (unsigned ilag)
  {
    if (ilag == 0)
      return get_covariance();
    
    return sample::get_crosscovariance (A, ilag, sample_size)
      + sample::get_crosscovariance (B, ilag, sample_size);
  }
};

/***************************************************************************
 *
 *  a superposition of two sources of electromagnetic radiation
 *
 ***************************************************************************/

class superposed : public combination
{
  Stokes<double> get_Stokes ();
  Vector<4, double> get_mean ();
  Matrix<4,4, double> get_covariance ();
};

/***************************************************************************
 *
 *  a composition of two sources of electromagnetic radiation
 *
 ***************************************************************************/

class composite : public combination
{
  double A_fraction;

public:

  composite (double fraction) { A_fraction = fraction; }

  Stokes<double> get_Stokes ();
  Vector<4, double> get_mean ();
  Matrix<4,4, double> get_covariance ();
};

/***************************************************************************
 *
 *  a disjoint combination of two sources of electromagnetic radiation
 *
 ***************************************************************************/

class disjoint : public combination
{
  double A_fraction;

public:

  disjoint (double fraction) { A_fraction = fraction; }

  Stokes<double> get_Stokes ();
  Vector<4, double> get_mean ();
  Matrix<4,4, double> get_covariance ();
  Matrix<4,4, double> get_crosscovariance (unsigned ilag);
};


/***************************************************************************
 *
 *  a post-detection boxcar-smoothed source of electromagnetic radiation
 *
 ***************************************************************************/

class boxcar_sample : public single
{
  std::vector< Stokes<double> > instances;
  unsigned smooth;
  unsigned current;

  void setup()
  {
    current = 0;
    instances.resize (smooth);
    for (unsigned i=1; i<smooth; i++)
      instances[i] = single::get_Stokes_instance();
  }

public:

  boxcar_sample (mode* s, unsigned n) : single(s) { smooth = n; }

  Stokes<double> get_Stokes_instance ()
  {
    if (instances.size() < smooth)
      setup ();

    instances[current] = single::get_Stokes_instance();
    current = (current + 1) % smooth;

    Stokes<double> result;
    for (unsigned i=0; i<smooth; i++)
      result += instances[i];

    result /= smooth;
    return result;
  }
};


class coherent : public combination
{
  Spinor<double> a;
  Spinor<double> b;

  field_transformer* a_xform;
  field_transformer* b_xform;

  mode* coupling;
  double coherence;
  
  bool built;
  void build();

 
public:

  coherent (double coherence);
  void set_normal (BoxMuller*);
  
  Stokes<double> get_Stokes ();
  Vector<4, double> get_mean ();
  Matrix<4,4, double> get_covariance ();
 
};


#endif
