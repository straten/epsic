/***************************************************************************
 *
 *   Copyright (C) 2024 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "true_math.h"
#include <math.h>

// see https://stackoverflow.com/questions/61941592/how-to-disable-fast-math-for-a-header-file-function

int true_signbit_float (float x)
{
  return signbit (x);
}

int true_signbit_double (double x)
{
  return signbit (x);
}

int true_signbit_long_double (long double x)
{
  return signbit (x);
}

int true_finite_float (float x)
{
  return isfinite (x) && !isnan (x);
}

int true_finite_double (double x)
{
  return isfinite (x) && !isnan (x);
}

int true_finite_long_double (long double x)
{
  return isfinite (x) && !isnan (x);
}
