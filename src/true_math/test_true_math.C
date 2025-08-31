/***************************************************************************
 *
 *   Copyright (C) 2024 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "true_math.h"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

int main ()
{
  const unsigned ntest = 3;

  /*
    Test true_math::finite
  */
  double some_nans[ntest];

  some_nans[0] = 1.0/0.0;
  some_nans[1] = 0.0/0.0;
  some_nans[2] = std::nan("");

  for (unsigned i = 0; i < ntest; i++ )
  {
    double x = some_nans[i];
    if (isfinite(x))
      cerr << "isfinite(" << x << ") fails" << endl;
    if (true_math::finite(x))
    {
      cerr << "true_math::finite(" << x << ") fails too!" << endl;
      return -1;
    }
    cerr << "true_math::finite(" << x << ") returns false as expected" << endl;
  }

  /*
    Test true_math::signbit
  */

  float zero = 0.0;
  float neg_zero = -0.0;

  if (true_math::signbit(zero))
  {
    cerr << "true_math::signbit(" << zero << ") fails" << endl;
    return -1;
  }
  cerr << "true_math::signbit(" << zero << ") returns false as expected" << endl;

  if (!true_math::signbit(neg_zero))
  {
    cerr << "true_math::signbit(" << neg_zero << ") fails" << endl;
    return -1;
  }
  cerr << "true_math::signbit(" << neg_zero << ") returns true as expected" << endl;

  return 0;
}

