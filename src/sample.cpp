
/***************************************************************************
 *
 *   Copyright (C) 2016 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "sample.h"

//! Sums the sample_size by sample_size square on the diagonal
/*! worker function for sub-classes */
Matrix<4,4, double> epsic::sample::get_covariance (mode* s, unsigned sample_size)
{
  Matrix<4,4, double> result = s->get_covariance ();

  // sum the sample_size instances along the diagonal
  result *= sample_size;

  for (unsigned ilag=1; ilag < sample_size; ilag++)
  {
    Matrix<4,4, double> out = s->get_crosscovariance (ilag);

#ifdef _DEBUG
    std::cerr << "ilag=" << ilag << " C=" << out << std::endl;
#endif

    /*
      multiply by two to also sum the symmetric lower triangle
      
      multiply by sample_size-ilag samples along the diagonal
      defined by i,i+ilag
    */
    out *= 2.0 * (sample_size-ilag);
    result += out;
  }
  
  /*          
    divide by sample_size squared because this function returns the
    covariance matrix of the sample mean
  */
    
  result /= sample_size * sample_size;
  return result;
}

//! Sums the sample_size by sample_size square off the diagonal
/*! starts at at_lag * sample_size off the diagonal */
Matrix<4,4, double> epsic::sample::get_crosscovariance (mode* s, unsigned at_lag,
						 unsigned sample_size)
{
  Matrix<4,4, double> result (0);
    
  for (unsigned ilag=0; ilag < sample_size; ilag++)
    for (unsigned jlag=0; jlag < sample_size; jlag++)
    {
      unsigned mode_lag = std::abs(int(at_lag*sample_size+ilag) - int(jlag));
      result += s->get_crosscovariance (mode_lag);
    }
  
  /*          
    divide by sample_size squared because this function returns the
    cross covariance matrix of the sample mean
  */
    
  result /= sample_size * sample_size;
  return result;
}
