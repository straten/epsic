/***************************************************************************
 *
 *   Copyright (C) 2020 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "covariant.h"

#include <assert.h>

using namespace std;

double covariant_mode::modulation ()
{
  if (amps.size() == 0)
    coordinator->get();

  double retval = amps.front();
  amps.pop();
  return retval;
}

covariant_coordinator::covariant_coordinator (double _correlation)
{
  correlation = _correlation;
  out[0] = out[1] = 0;
}

modulated_mode* covariant_coordinator::get_modulated_mode (unsigned index, mode* in)
{
  assert (index < 2);

  if (!out[index])
  {
    if (!in)
      throw std::runtime_error( "covariant_coordinator::get_modulated_mode "
                                "input not set" );

    out[index] = new covariant_mode ( in );
    out[index]->coordinator = this;
  }

  return out[index];
}

void covariant_coordinator::get()
{
  for (unsigned i=0; i<2; i++)
    if (!out[i])
      throw std::runtime_error( "covariant_coordinator::get "
                                "output not set" );

  double a0, a1;
  get_modulation (a0, a1);

  // cerr << "covariant_coordinator::get a0=" << a0 << " a1=" << a1 << endl;

  out[0]->amps.push( a0 );
  out[1]->amps.push( a1 );
}

Matrix<2,2,double> sqrt (const Matrix<2,2,double>& C)
{
  double det = C[0][0]*C[1][1] - C[0][1]*C[1][0];
  double trace = C[0][0]+C[1][1];

  double s = sqrt(det);
  double t = sqrt(trace + 2*s);

  Matrix<2,2,double> result (s);
  result += C;
  result /= t;
  return result;
}

bivariate_lognormal_modes::~bivariate_lognormal_modes ()
{
  mean /= count;
  meansq /= count;
  meansq -= outer (mean, mean);

  cerr << "\n"
          "bivariate_lognormal_modes mean=" << mean << " "
          "rho=" << meansq[0][1]/sqrt(meansq[0][0]*meansq[1][1]) << "\n"
          "covar=\n" << meansq << endl;
}

void bivariate_lognormal_modes::build ()
{
  double correlation = get_correlation();

  Matrix<2,2,double> covar;

  covar[0][0] = log_sigma[0]*log_sigma[0];
  covar[1][1] = log_sigma[1]*log_sigma[1];

  double beta0 = sqrt( exp(covar[0][0]) - 1.0 );
  double beta1 = sqrt( exp(covar[1][1]) - 1.0 );

  double denom = beta0 * beta1;
  double max_correlation = (exp(log_sigma[0]*log_sigma[1]) - 1.0) / denom;
  double min_correlation = (exp(-log_sigma[0]*log_sigma[1]) - 1.0) / denom;

  if (correlation > max_correlation)
  {
    cerr << "bivariate_lognormal_modes::build correlation=" << correlation
         << " > max=" << max_correlation << endl;

    throw std::runtime_error( "bivariate_lognormal_modes::build "
                              "maximum correlation exceeded" );
  }

  if (correlation < min_correlation)
  { 
    cerr << "bivariate_lognormal_modes::build correlation=" << correlation
         << " < min =" << min_correlation << endl;
    
    throw std::runtime_error( "bivariate_lognormal_modes::build "
                              "minimum correlation exceeded" );
  }

  covar[0][1] = covar[1][0] = log( correlation * beta0 * beta1 + 1 );

  correlator = sqrt(covar);

  meansq = 0;
  mean = 0;
  count = 0;

  built = true;
}

void bivariate_lognormal_modes::get_modulation (double& A, double& B)
{
  if (!built)
    build ();

  assert (normal != NULL);

  Vector<2,double> amps;
  amps[0] = normal->evaluate();
  amps[1] = normal->evaluate();

  amps = correlator * amps;

  amps[0] = A = exp (amps[0] - 0.5*log_sigma[0]*log_sigma[0]);
  amps[1] = B = exp (amps[1] - 0.5*log_sigma[1]*log_sigma[1]);

  mean += amps;
  meansq += outer (amps, amps);
  count ++;
}

void bivariate_lognormal_modes::set_beta (unsigned index, double beta)
{
  assert (index < 2);
  log_sigma[index] = sqrt( log( beta*beta + 1.0 ) );
  built = false;
}

