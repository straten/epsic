//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2016 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/modulated.h

#ifndef __modulated_H
#define __modulated_H

/***************************************************************************
 *
 *  an amplitude modulated mode of electromagnetic radiation
 *
 ***************************************************************************/

class modulated_mode : public mode_decorator
{
#if _DEBUG
  double tot, totsq;
  uint64_t count;
#endif

public:

  modulated_mode (mode* s) : mode_decorator(s) { 
#if -DEBUG
    tot=0; totsq=0; count=0; 
#endif
}

  // return a random scalar modulation factor
  virtual double modulation () = 0;

  // return the mean of the scalar modulation factor
  virtual double get_mod_mean () = 0;

  // return the variance of the scalar modulation factor
  virtual double get_mod_variance () = 0;

  Spinor<double> get_field ()
  {
    double mod = modulation();
#if _DEBUG
    tot+=mod;
    totsq+=mod*mod;
    count+=1;
#endif
    return sqrt(mod) * source->get_field();
  }

  Matrix<4,4,double> get_covariance ()
  {
    double mean = get_mod_mean();
    double var = get_mod_variance();

#if _DEBUG
    cerr << "modulated_mode::get_covariance expected mean=" << mean
	 << " var=" << var << endl;
    tot /= count;
    totsq /= count;
    totsq -= tot*tot;
    cerr << " measured mean=" << tot << " var=" << totsq << endl;
#endif

    Matrix<4,4,double> C = source->get_covariance();
    C *= (mean*mean + var);
    Matrix<4,4,double> o = outer (get_Stokes(), get_Stokes());
    o *= var;
    return C + o;
  }

  Stokes<double> get_mean ()
  {
    return get_mod_mean () * source->get_mean();
  }

};

class lognormal_mode : public modulated_mode
{
  // standard deviation of the logarithm of the random variate
  double log_sigma;

public:

  lognormal_mode (mode* s, double ls) : modulated_mode (s) { log_sigma = ls;}

  // return a random scalar modulation factor
  double modulation ()
  {
    return exp ( log_sigma * (get_normal()->evaluate() - 0.5*log_sigma) ) ;
  }

  double get_mod_mean () { return 1.0; }
  
  double get_mod_variance () { return exp(log_sigma*log_sigma) - 1.0; }

};


class boxcar_modulated_mode : public modulated_mode
{
  std::vector< double > instances;
  unsigned smooth;
  unsigned current;

  void setup()
  {
    current = 0;
    instances.resize (smooth);
    for (unsigned i=1; i<smooth; i++)
      instances[i] = mod->modulation();
  }

  modulated_mode* mod;

public:

  boxcar_modulated_mode (modulated_mode* s, unsigned n)
    : modulated_mode(s->get_source()) { smooth = n; mod = s; }

  double modulation ()
  {
    if (instances.size() < smooth)
      setup ();

    instances[current] = mod->modulation();
    current = (current + 1) % smooth;

    double result = 0.0;
    for (unsigned i=0; i<smooth; i++)
      result += instances[i];

    result /= smooth;

    return result;
  }

  double get_mod_variance ()
  {
    return mod->get_mod_variance() / smooth;
  }

  double get_mod_mean ()
  {
    return mod->get_mod_mean ();
  }

  //! Return the sum of the intensity autocorrelation function
  double get_autocorrelation (unsigned nsample) 
  {
    double nsum = std::min (smooth, nsample);
    double result = 0;

    // sum all of the elements in the upper triangle of the covariance matrix
    for (unsigned offset=1; offset < nsum; offset++)
    {
      double acf = (smooth - offset) / double (smooth * smooth);
      result += (nsample - offset) * acf;
    }

    // multiply by two to also sum the symmetric lower triangle
    return 2.0 * result * mod->get_mod_variance();
  }
};


class square_modulated_mode : public modulated_mode
{
  unsigned width;
  unsigned current;
  double value;

  modulated_mode* mod;

public:

  square_modulated_mode (modulated_mode* s, unsigned n)
    : modulated_mode(s->get_source()) { width = n; current = n; mod = s; }

  double modulation ()
  {
    if (current == width)
    {
      value = mod->modulation();
      current = 0;
    }

    current ++;
    return value;
  }

  double get_mod_variance ()
  {
    return mod->get_mod_variance();
  }

  double get_mod_mean ()
  {
    return mod->get_mod_mean ();
  }

  //! Return the sum of the intensity autocorrelation function
  double get_autocorrelation (unsigned nsample) 
  {
    // sum the width X width squares along the diagonal
    double result = 0;
    while (nsample > width)
    {
      result += (width-1)*width;
      nsample -= width;
    }

    if (nsample)
      result += (nsample-1)*nsample;
    
    return result * mod->get_mod_variance();
  }
};

#endif
