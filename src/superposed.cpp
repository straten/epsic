/***************************************************************************
 *
 *   Copyright (C) 2016 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "sample.h"

Stokes<double> epsic::superposed::get_Stokes ()
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

Vector<4, double> epsic::superposed::get_mean ()
{
  return A->get_mean() + B->get_mean();
}

//! Implements Equation (42) of van Straten & Tiburzi (2017)
Matrix<4,4, double> epsic::superposed::get_covariance ()
{
  Matrix<4,4, double> result = sample::get_covariance (A, sample_size);
  result += sample::get_covariance (B, sample_size);
  
  Stokes<double> mean_A = A->get_mean();
  Stokes<double> mean_B = B->get_mean();
  
  /* Minkowski::outer implements A \otimes B - 0.5 \eta A \cdot B
     such that Minkowski::outer(A,B) +  Minkowski::outer(B,A)
     yields Equation (43) of van Straten & Tiburzi (2017) */
  Matrix<4,4, double> xcovar = Minkowski::outer(mean_A, mean_B);

  // assumes that the mean of the modulating function is unity
  xcovar *= (1.0 + intensity_covariance) / sample_size;

  Matrix<4,4, double> extra = outer(mean_A, mean_B);
  extra *= intensity_covariance / sample_size;

  xcovar += extra;

  result += xcovar + transpose(xcovar);

  return result;
}

