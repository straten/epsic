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

  // worker function for sub-classes
  Matrix<4,4, double> get_covariance (mode* s, unsigned n)
  {
    std::cerr << "sample::get_covariance n=" << n << std::endl;

    Matrix<4,4, double> result = s->get_covariance ();

    for (unsigned ilag=1; ilag < n; ilag++)
    {
      Matrix<4,4, double> out = s->get_crosscovariance (ilag);

      std::cerr << "ilag=" << ilag << " C=" << out << std::endl;

      // multiply by two to also sum the symmetric lower triangle
      out *= 2.0*double(n-ilag)/double(n);
      result += out;
    }
    
    result /= n;
    return result;
  }
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
};


/***************************************************************************
 *
 *  a combination of two sources of electromagnetic radiation
 *
 ***************************************************************************/

class combination : public sample
{
public:
  mode* A;
  mode* B;

  combination () { A = new mode; B = new mode; }

  void set_normal (BoxMuller* n) { A->set_normal(n); B->set_normal(n); }
};

/***************************************************************************
 *
 *  a superposition of two sources of electromagnetic radiation
 *
 ***************************************************************************/

class superposed : public combination
{
  Stokes<double> get_Stokes ()
  {
    Stokes<double> result;
    for (unsigned i=0; i<sample_size; i++)
    {
      Spinor<double> e_A = A->get_field();
      Spinor<double> e_B = B->get_field();

      Vector<4, double> tmp;
      compute_stokes (tmp, e_A + e_B);
      result += tmp;
    }
    result /= sample_size;
    return result;
  }

  Vector<4, double> get_mean ()
  {
    return A->get_mean() + B->get_mean();
  }

  Matrix<4,4, double> get_covariance ()
  {
    Matrix<4,4, double> result = sample::get_covariance (A, sample_size);
    result += sample::get_covariance (B, sample_size);

    Stokes<double> mean_A = A->get_mean();
    Stokes<double> mean_B = B->get_mean();
    Matrix<4,4, double> xcovar = Minkowski::outer(mean_A, mean_B);
    xcovar /= sample_size;
    result += xcovar + transpose(xcovar);
    
    return result;
  }
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

  Stokes<double> get_Stokes ()
  {
    unsigned A_sample_size = A_fraction * sample_size;

    Stokes<double> result;

    for (unsigned i=0; i<sample_size; i++)
    {
      Spinor<double> e;
      if (i < A_sample_size)
	e = A->get_field();
      else
	e = B->get_field();

      Vector<4, double> tmp;
      compute_stokes (tmp, e);
      result += tmp;
    }

    result /= sample_size;
    return result;
  }

  Vector<4, double> get_mean ()
  {
    unsigned A_sample_size = A_fraction * sample_size;
    unsigned B_sample_size = sample_size - A_sample_size;
    Vector<4,double> result 
      = A_sample_size * A->get_mean() + B_sample_size * B->get_mean();
    result /= sample_size;
    return result;
  }

  Matrix<4,4, double> get_covariance ()
  {
    unsigned A_sample_size = A_fraction * sample_size;
    unsigned B_sample_size = sample_size - A_sample_size;

    Matrix<4,4,double> C_A = sample::get_covariance (A, A_sample_size);
    Matrix<4,4,double> C_B = sample::get_covariance (B, B_sample_size);

    // A_fraction * sample_size may not be an integer number of instances
    double f_A = A_sample_size / double(sample_size);
    double f_B = B_sample_size / double(sample_size);

    C_A *= f_A * f_A;
    C_B *= f_B * f_B;

    return C_A + C_B;
  }
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

  Stokes<double> get_Stokes ()
  {
    bool mode_A = random_double() < A_fraction;
    mode* e = (mode_A) ? A : B;

    Stokes<double> result;

    for (unsigned i=0; i<sample_size; i++)
    {
      Vector<4, double> tmp;
      compute_stokes (tmp, e->get_field());
      result += tmp;
    }

    result /= sample_size;
    return result;
  }

  Vector<4, double> get_mean ()
  {
    return A_fraction * A->get_mean() + (1-A_fraction) * B->get_mean();
  }

  Matrix<4,4, double> get_covariance ()
  {
    Matrix<4,4,double> C_A = sample::get_covariance (A, sample_size);
    Matrix<4,4,double> C_B = sample::get_covariance (B, sample_size);

    Vector<4,double> diff = A->get_mean() - B->get_mean();
    Matrix<4,4,double> D = outer (diff, diff);

    C_A *= A_fraction;
    C_B *= (1-A_fraction);
    D *= A_fraction * (1-A_fraction);

    return C_A + C_B + D;
  }
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

#endif
