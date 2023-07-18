#pragma GCC optimize ("no-fast-math")
#include <math.h>

// see https://stackoverflow.com/questions/61941592/how-to-disable-fast-math-for-a-header-file-function

int myfinite (double x)
{
  return isfinite (x) && !isnan (x);
}

