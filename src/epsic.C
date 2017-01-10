/***************************************************************************
 *
 *   Copyright (C) 2016 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

/*

  Simulates the polarization of electromagnetic radiation.

*/

#include "Pauli.h"
#include "Dirac.h"
#include "Jacobi.h"

#include "mode.h"
#include "modulated.h"
#include "smoothed.h"
#include "sample.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <inttypes.h>
#include <exception>

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
  while ((c = getopt(argc, argv, "b:dr:hn:N:m:M:s:SC:D:A:B:l:opRX")) != -1)
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

    case 'b':
      if (optarg[0]=='B')
	setup_B.smooth_modulator = atoi (optarg+1);
      else
	setup_A.smooth_modulator = atoi (optarg);
      break;

    case 'r':
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

    case 'R':
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

