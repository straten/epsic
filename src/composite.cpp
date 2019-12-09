/***************************************************************************
 *
 *   Copyright (C) 2016 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "sample.h"

Stokes<double> composite::get_Stokes ()
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

Vector<4, double> composite::get_mean ()
{
  unsigned A_sample_size = A_fraction * sample_size;
  unsigned B_sample_size = sample_size - A_sample_size;
  Vector<4,double> result =
    A_sample_size * A->get_mean() +
    B_sample_size * B->get_mean();
  result /= sample_size;
  return result;
}

//! Implements Equation (59) of van Straten & Tiburzi (2017)
Matrix<4,4, double> composite::get_covariance ()
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

