
#include "myfinite.h"

#include <iostream>
#include <vector>
#include <cmath>

// see https://stackoverflow.com/questions/61941592/how-to-disable-fast-math-for-a-header-file-function

using namespace std;

int main ()
{
  const unsigned ntest = 3;

  double some_nans[ntest];

  some_nans[0] = 1.0/0.0;
  some_nans[1] = 0.0/0.0;
  some_nans[2] = std::nan("");

  for (unsigned i = 0; i < ntest; i++ )
  {
    double x = some_nans[i];
    if (isfinite(x))
      cerr << "isfinite(" << x << ") fails" << endl;
    if (myfinite(x))
    {
      cerr << "myfinite(" << x << ") fails too!" << endl;
      return -1;
    }
    cerr << "myfinite(" << x << ") returns false as expected" << endl;
  }

  return 0;
}

