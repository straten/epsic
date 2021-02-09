/***************************************************************************
 *
 *   Copyright (C) 2003 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
#include "Vector.h"
#include "random.h"

using namespace std;

int main ()
{
  Vector<3,double> v1;
  Vector<3,double> v2;

  random_vector (v1, 10.0);
  random_vector (v2, 10.0);

  cerr << "v1=" << v1 << " v2=" << v2 << " v1+v2=" << v1+v2 
       << " v1*v2=" << v1*v2 << endl;

  return 0;
}

