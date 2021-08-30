/***************************************************************************
 *
 *   Copyright (C) 2016 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "sample.h"

using namespace std;

void add (Stokes<double>& result, Spinor<double>& e)
{
  Vector<4, double> tmp;
  compute_stokes (tmp, e);
  result += tmp;
}

Stokes<double> composite::get_Stokes ()
{
  unsigned A_sample_size = A_fraction * sample_size;
  unsigned B_sample_size = sample_size - A_fraction;
  unsigned max_size = std::max (A_sample_size, B_sample_size);
  
  Stokes<double> result;
  Spinor<double> e;

  // ensure that A->get_field and B->get_field are called equal number of times 
  for (unsigned i=0; i<max_size; i++)
  {
    e = A->get_field();
    if (i < A_sample_size)
      add (result, e);

    e = B->get_field();
    if (i < B_sample_size)
      add (result, e);
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

  Stokes<double> mean_A = A->get_mean();
  Stokes<double> mean_B = B->get_mean();

  Matrix<4,4, double> extra = outer(mean_A, mean_B);
  extra *= std::min(f_A,f_B) * intensity_covariance / sample_size;

  return C_A + C_B + extra + transpose(extra);
}

