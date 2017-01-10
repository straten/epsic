
#include "random.h"

#include <stdlib.h>
#include <sys/time.h>

// #define _DEBUG 1

#if _DEBUG
#include <iostream>
using namespace std;
#endif

static const double random_max = RAND_MAX;

uint64_t usec_seed ()
{
  struct timeval t;
  gettimeofday (&t, NULL);
  // cerr << "usec_seed: usec=" << t.tv_usec << endl;
  return t.tv_usec;
}

void random_init ()
{
  srandom (usec_seed());
}

double random_double ()
{
  double val = double(random()) / random_max;
#if _DEBUG
  static double vmin = 1000;
  vmin = std::min(vmin,val);
  static double vmax = -1000;
  vmax = std::max(vmax,val);
  cerr << "check " << vmin << " " << vmax << endl;
#endif
  return val;
}
 
