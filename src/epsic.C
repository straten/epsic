/***************************************************************************
 *
 *   Copyright (C) 2008 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

/*

  Simulate polarized noise

*/

#include "Minkowski.h"
#include "Stokes.h"
#include "Jones.h"
#include "Pauli.h"
#include "Dirac.h"
#include "Jacobi.h"

#include "BoxMuller.h"
#include "random.h"

#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <inttypes.h>

// #define _DEBUG 1

using namespace std;

void usage ()
{
  cout <<

    "simpol: simulate polarized noise and compute statistics \n" 
    " \n"
    "options: \n"
    " \n"
    " -N Msamp    simulate Msamp mega samples \n"
    " -n Nint     integrate Nint samples before further processing \n"
    " -m Nsamp    box-car smooth over Nsamp samples before detection \n"
    " -M Nsamp    box-car smooth over Nsamp samples after detection \n"
    " -s i,q,u,v  single source with specified Stokes parameters \n"
    " -S          superposed modes \n"
    " -X          cross-correlate the modes after detection \n"
    " -C f_A      composite modes with fraction of instances in mode A \n"
    " -D F_A      disjoint modes with fraction of samples in mode A \n"
    " -A i,q,u,v  set the Stokes parameters of mode A \n"
    " -B i,q,u,v  set the Stokes parameters of mode B \n"
    " -l sigma    modulate mode A using a log-normal variate \n"
    " -a Nsamp    box-car smooth the amplitude modulation function \n"
    " -g Nsamp    use square impulse amplitude modulation function \n"
    " -p          print real part of x and y, plus diff. phase \n"
    " -r          print the statistics of the coherency matrix \n"
    " -d          print only the variances of each Stokes parameter \n"
       << endl;
}

/***************************************************************************
 *
 *  spinor class represents an instance of the electric field
 *
 ***************************************************************************/

template<typename T>
class spinor
{
public:
  spinor (const std::complex<T>& _x, const std::complex<T>& _y) : x(_x), y(_y){}
  spinor () { }

  std::complex<T> x;
  std::complex<T> y;

  const spinor& operator *= (T scale) { x *= scale; y *= scale; return *this; }
  const spinor& operator /= (T norm) { x /= norm; y /= norm; return *this; }
  const spinor& operator += (const spinor& e) { x+=e.x; y+=e.y; return *this; }
};

template <typename T>
ostream& operator << (ostream& os, const spinor<T>& s)
{
  os << s.x << " " << s.y;
  return os;
}

template <typename T>
const spinor<T> operator + (spinor<T> s, const spinor<T>& t)
{
  return s += t;
}

//! returns the input spinor transformed by the Jones matrix
template<typename T>
const spinor<T> operator * (const Jones<T>& j, const spinor<T>& in)
{
  return spinor<T> ( j.j00 * in.x + j.j01 * in.y,
		     j.j10 * in.x + j.j11 * in.y );
}

template<typename T>
const spinor<T> operator * (double a, spinor<T> in)
{
  return in *= a;
}

template<typename T>
const spinor<T> operator * (spinor<T> in, double a)
{
  return in *= a;
}

/***************************************************************************
 *
 *  a single source of electromagnetic radiation
 *
 ***************************************************************************/

class mode
{
public:
  mode ();
  virtual ~mode () { }

  virtual void set_Stokes (const Stokes<double>& mean);
  virtual Stokes<double> get_Stokes () { return mean; }

  //! Return the expected mean Stokes parameters
  virtual Stokes<double> get_mean () { return mean; }

  //! Return the expected covariances between the Stokes parameters
  virtual Matrix<4,4, double> get_covariance ()
  { return Minkowski::outer (mean, mean); }

  //! Return the sum of the intensity autocorrelation function
  virtual double get_autocorrelation (unsigned nsample) { return 0; }

  //! Return a random instance of the electric field vector
  virtual spinor<double> get_field ();

  //! Return BoxMuller object used to generate normally distributed numbers
  virtual BoxMuller* get_normal () { return normal; }
  virtual void set_normal (BoxMuller* n) { normal = n; }

  void set_power_law (double _inner_scale) { inner_scale = _inner_scale; }
  void powerlaw (spinor<double>&);

private:
  Stokes<double> mean;
  Jones<double> polarizer;

  BoxMuller* normal;
  double rms;

  double inner_scale;
};
  
mode::mode ()
{
  normal = 0;
  rms = 0.5;

  inner_scale = 0.0;

  set_Stokes (Stokes<double>(1.0));
}

void mode::set_Stokes (const Stokes<double>& _mean)
{
  mean = _mean;

  Quaternion<double,Hermitian> root = sqrt (natural(mean));
  polarizer = convert (root);
}

spinor<double> mode::get_field ()
{
  if (!normal)
    throw Error (InvalidState, "mode::get_field", "BoxMuller not set");

  BoxMuller& gasdev = *normal;

  complex<double> x (rms * gasdev(), rms * gasdev());
  complex<double> y (rms * gasdev(), rms * gasdev());

  spinor<double> e (x, y);

  if (inner_scale)
    powerlaw (e);

  return polarizer * e;
}

void mode::powerlaw (spinor<double>& e)
{
  double y = random_double(); 
  double index = -4;

  e *= pow( pow(inner_scale,index+1) * (1.0 - y), 1.0/(index+1) );
}

template<typename T, typename U>
void compute_stokes (Vector<4,T>& stokes, const spinor<U>& e)
{
  double var_x = norm(e.x);
  double var_y = norm(e.y);

  complex<double> c_xy = conj(e.x) * e.y;
  
  stokes[0] = var_x + var_y;
  stokes[1] = var_x - var_y;
  stokes[2] = 2.0*c_xy.real();
  stokes[3] = 2.0*c_xy.imag();
}

class mode_decorator : public mode
{
protected:
  mode* source;

public:
  mode_decorator (mode* s) { source = s; }
  mode* get_source () { return source; }

  void set_Stokes (const Stokes<double>& mean) { source->set_Stokes(mean); }
  Stokes<double> get_Stokes () { return source->get_Stokes(); }

  Matrix<4,4, double> get_covariance () { return source->get_covariance(); }
  Stokes<double> get_mean () { return source->get_mean(); }

  spinor<double> get_field () { return source->get_field(); }
  BoxMuller* get_normal () { return source->get_normal(); }
  void set_normal (BoxMuller* n) { source->set_normal(n); }
};

/***************************************************************************
 *
 *  a boxcar-smoothed mode of electromagnetic radiation
 *
 ***************************************************************************/

class boxcar_mode : public mode_decorator
{
  vector< spinor<double> > instances;
  unsigned smooth;
  unsigned current;

  void setup()
  {
    current = 0;
    instances.resize (smooth);
    for (unsigned i=1; i<smooth; i++)
      instances[i] = source->get_field();
  }

public:

  boxcar_mode (mode* s, unsigned n) : mode_decorator(s) { smooth = n; }

  spinor<double> get_field ()
  {
    if (instances.size() < smooth)
      setup ();

    instances[current] = source->get_field();
    current = (current + 1) % smooth;

    spinor<double> result;
    for (unsigned i=0; i<smooth; i++)
      result += instances[i];

    result /= sqrt(smooth);
    return result;
  }
};

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

  spinor<double> get_field ()
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
  vector< double > instances;
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
    Matrix<4,4, double> result = s->get_covariance ();

    double acf = s->get_autocorrelation (n);
    Matrix<4,4, double> out = outer( s->get_mean(), s->get_mean() );
    out *= acf / n;
    result += out;
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
    spinor<double> e = source->get_field();
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
      spinor<double> e_A = A->get_field();
      spinor<double> e_B = B->get_field();

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
      spinor<double> e;
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
  vector< Stokes<double> > instances;
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

class mode_setup
{
public:
  // box-car smoothing width pre-detection
  unsigned smooth_before;
  // box-car smoothing of modulation function
  unsigned smooth_modulator;
  // width of square modulation function
  unsigned square_modulator;
  // variance of logarithm of modulation function
  double log_sigma;

  mode_setup ()
  {
    smooth_before = 0;
    smooth_modulator = 0;
    square_modulator = 0;
    log_sigma = 0;
  }

  mode* setup_mode (mode* s)
  {
    modulated_mode* mod = 0;

    if (log_sigma)
      s = mod = new lognormal_mode (s, log_sigma);

    if (smooth_modulator > 1 && mod)
      s = new boxcar_modulated_mode (mod, smooth_modulator);

    if (square_modulator > 1 && mod)
      s = new square_modulated_mode (mod, square_modulator);

    if (smooth_before > 1)
      s = new boxcar_mode (s, smooth_before);

    return s;
  }
};

double sqr (double x) { return x*x; }

int main (int argc, char** argv)
{
  uint64_t ndat = 1024 * 1024;     // number of Stokes samples
  unsigned nint = 1;               // number of instances in Stokes sample

  unsigned smooth_after = 0;       // box-car smoothing width post-detection

  bool verbose = false;

  Stokes<double> stokes = 1.0;
  bool subtract_outer_population_mean = false;

  mode source;
  combination* dual = NULL;
  sample* stokes_sample = NULL;

  mode_setup setup_A;
  mode_setup setup_B;

  bool cross_correlate = false;

  bool print = false;
  bool rho_stats = false;
  bool variances_only = false;

  int c;
  while ((c = getopt(argc, argv, "a:dg:hn:N:m:M:s:SC:D:A:B:l:oprX")) != -1)
  {
    switch (c)
    {

    case 'h':
      usage ();
      return 0;

    case 's':
    case 'A':
    case 'B':
    {
      double i,q,u,v;
      if (sscanf (optarg, "%lf,%lf,%lf,%lf", &i,&q,&u,&v) != 4)
      {
	cerr << "Error parsing " << optarg << " as Stokes 4-vector" << endl;
	return -1;
      }
      stokes = Stokes<double> (i,q,u,v);

      if (stokes.abs_vect() > i) {
	cerr << "Invalid Stokes parameters (p>I) " << stokes << endl;
	return -1;
      }

      if (dual && c=='A')
	dual->A->set_Stokes( stokes );
      if (dual && c=='B')
	dual->B->set_Stokes( stokes );

      break;
    }

    case 'l':
      if (optarg[0]=='B')
	setup_B.log_sigma = atof (optarg+1);
      else
	setup_A.log_sigma = atof (optarg);
      break;
      
    case 'N':
      ndat = ndat * atof (optarg);
      break;

    case 'n':
      nint = atoi (optarg);
      break;

    case 'M':
      smooth_after = atoi (optarg);
      break;

    case 'm':
      if (optarg[0]=='B')
	setup_B.smooth_before = atoi (optarg+1);
      else
	setup_A.smooth_before = atoi (optarg);
      break;

    case 'a':
      if (optarg[0]=='B')
	setup_B.smooth_modulator = atoi (optarg+1);
      else
	setup_A.smooth_modulator = atoi (optarg);
      break;

    case 'g':
      if (optarg[0]=='B')
	setup_B.square_modulator = atoi (optarg+1);
      else
	setup_A.square_modulator = atoi (optarg);
      break;

    case 'o':
      subtract_outer_population_mean = true;
      break;

    case 'p':
      print = true;
      break;

    case 'S':
      dual = new superposed;
      break;

    case 'X':
      dual = new superposed;
      cross_correlate = true;
      break;

    case 'C':
      dual = new composite( atof(optarg) );
      break;

    case 'D':
      dual = new disjoint( atof(optarg) );
      break;

    case 'r':
      rho_stats = true;
      break;

    case 'd':
      variances_only = true;
      break;

    case 'v':
      verbose = true;
      break;
    }
  }

  cerr << "Simulating " << ndat << " Stokes samples" << endl;

  random_init ();
  BoxMuller gasdev (time(NULL));

  uint64_t ntot = 0;
  double totp = 0;
  Vector<4, double> tot;
  Matrix<4,4, double> totsq;

  Matrix<2,2, complex<double> > tot_rho;
  Matrix<4,4, complex<double> > totsq_rho;

  source.set_Stokes (stokes);

  if (dual)
  {
    dual->set_normal (&gasdev);
    stokes_sample = dual;

    dual->A = setup_A.setup_mode (dual->A);
    dual->B = setup_B.setup_mode (dual->B);
  }
  else
  {
    source.set_normal (&gasdev);
    mode* s = setup_A.setup_mode(&source);

    if (smooth_after > 1)
      stokes_sample = new boxcar_sample (s, smooth_after);
    else
      stokes_sample = new single(s);
  }

  stokes_sample->sample_size = nint;

  for (uint64_t idat=0; idat<ndat; idat++)
  {
    Vector<4, double> mean_stokes;

    mean_stokes = stokes_sample->get_Stokes();

    tot += mean_stokes;
    totsq += outer(mean_stokes, mean_stokes);

    double psq = sqr(mean_stokes[1])+sqr(mean_stokes[2])+sqr(mean_stokes[3]);
    totp += sqrt(psq)/mean_stokes[0];
      
    if (rho_stats)
    {
      Matrix<2,2, complex<double> > rho = convert (Stokes<double>(mean_stokes));
    
      tot_rho += rho;
      totsq_rho += direct (rho, rho);
    }

    ntot ++;
  }

  totp /= ntot;
  tot /= ntot;
  totsq /= ntot;

  if (subtract_outer_population_mean)
    totsq -= outer(stokes,stokes);
  else
    totsq -= outer(tot,tot);

  if (variances_only)
  {
    for (unsigned i=0; i<4; i++)
      cout << "var[" << i << "] = " << totsq[i][i] << endl;
    
    return 0;
  }

  Vector<4, double> expected_mean;
  Matrix<4,4, double> expected_covariance;

  expected_mean = stokes_sample->get_mean ();
  expected_covariance = stokes_sample->get_covariance();

  cerr << "\n"
    " ******************************************************************* \n"
    "\n"
    " STOKES PARAMETERS \n"
    "\n"
    " ******************************************************************* \n"
       << endl;

  cerr << "mean sample dop=" << totp << endl << endl;

  cerr << "modulation index=" << sqrt(totsq[0][0])/tot[0] << endl << endl;

  cerr << "mean=" << tot << endl;
  cerr << "expected=" << expected_mean << endl;

  cerr << "\ncovar=\n" << totsq << endl;
  cerr << "expected=\n" << expected_covariance << endl;

  if (!rho_stats)
    return 0;

  cerr << "\n"
    " ******************************************************************* \n"
    "\n"
    " COHERENCY MATRIX \n"
    "\n"
    " ******************************************************************* \n"
       << endl;

  tot_rho /= ndat;
  totsq_rho /= ndat;

  cerr << "rho sq=\n" << totsq_rho << endl;

  totsq_rho -= direct(tot_rho,tot_rho);

  cerr << "rho mean=\n" << tot_rho << endl;
  cerr << "rho covar=\n" << totsq_rho << endl;

  Matrix<4,4, complex<double> > candidate;
  for (unsigned i=0; i<4; i++)
    for (unsigned j=0; j<4; j++)
      {
	Matrix<4,4,complex<double> > temp = Dirac::matrix (i,j);
	temp *= expected_covariance[i][j] * 0.25;

	candidate += temp;
      }

  cerr << "candidate=\n" << candidate << endl;

  Matrix<4, 4, complex<double> > eigenvectors;
  Vector<4, double> eigenvalues;

  Matrix<4, 4, complex<double> > temp = candidate;
  Jacobi (temp, eigenvectors, eigenvalues);

  for (unsigned i=0; i<4; i++)
    cerr << "e_" << i << "=" << eigenvalues[i] << "  v=" << eigenvectors[i] << endl;

  return 0;
}

