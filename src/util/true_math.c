/***************************************************************************
 *
 *   Copyright (C) 2024 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#pragma GCC optimize ("no-fast-math")
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

int true_finite_float (float x)
{
  return isfinite (x) && !isnan (x);
}

int true_finite_double (double x)
{
  return isfinite (x) && !isnan (x);
}

