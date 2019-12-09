
/***************************************************************************
 *
 *   Copyright (C) 2016 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "sample.h"

Stokes<double> disjoint::get_Stokes ()
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

Vector<4, double> disjoint::get_mean ()
{
  return A_fraction * A->get_mean() + (1-A_fraction) * B->get_mean();
}

//! Implements Equation (39) of van Straten & Tiburzi (2017)
Matrix<4,4, double> disjoint::get_covariance ()
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

Matrix<4,4, double> disjoint::get_crosscovariance (unsigned ilag)
{
  if (ilag == 0)
    return get_covariance();
  
  Matrix<4,4,double> Acov = A->get_crosscovariance(ilag);
  Acov *= A_fraction * A_fraction;
  
  Matrix<4,4,double> Bcov = B->get_crosscovariance(ilag);
  Bcov *= (1-A_fraction) * (1-A_fraction);
  
  return Acov + Bcov;
}
