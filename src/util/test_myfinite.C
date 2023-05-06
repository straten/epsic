
#include "myfinite.h"

#include <iostream>
#include <cmath>

// see https://stackoverflow.com/questions/61941592/how-to-disable-fast-math-for-a-header-file-function

using namespace std;

int main ()
{
  float x = 1.0 / 0.0;
  if (myfinite(x))
    cerr << "myfinite(1/0) fails" << endl;
  if (myfinite(x))
  {
    cerr << "myfinite(1/0) fails too!" << endl;
    return -1;
  }
  cerr << "myfinite(1/0) returns false as expected" << endl;
  return 0;
}

