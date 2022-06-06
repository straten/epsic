/***************************************************************************
 *
 *   Copyright (C) 2022 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

/* This test driver verifies that the rotation matrices returned by the
   rotation function rotate vectors about an axis in a manner that is
   consistent with the right-hand rule. */

#include "Matrix.h"

#include <iostream>
using namespace std;

static const char* name = "xyz";

void test_Rotation (double theta, double phi, unsigned axis, unsigned perm)
{
  Matrix<3,3,double> rot = rotation (Vector<3,double>::basis(axis), phi);

  unsigned i = (axis+perm+1)%3;
  unsigned j = (axis+perm+2)%3;
  unsigned k = (axis+perm)%3;

  Vector<3,double> vec;

  vec[i] = cos(theta);
  vec[j] = sin(theta);
  vec[k] = 0.2;

  vec = rot * vec;

  double result = theta + phi;

  if (fabs (vec[i] - cos(result)) > 1e-12)
  {
    cerr << name[i] << "=" << vec[i] << " != " << cos(result) << endl;
    exit(-1);
  }

  if (fabs (vec[j] - sin(result)) > 1e-12)
  {
    cerr << name[j] << "=" << vec[j] << " != " << sin(result) << endl;
    exit(-1);
  }

  if (fabs (vec[k] - 0.2) > 1e-12) {
    cerr << name[k] << "=" << vec[k] << " != " << 0.2 << endl;
    exit(-1);
  }

}
 
void test_loop (unsigned axis, unsigned perm)
{
  double increment = M_PI/87.0;

  for (double theta = -2*M_PI; theta < 2*M_PI+increment; theta += increment)
    for (double phi = -M_PI; phi < M_PI+increment; phi += increment)
      test_Rotation (theta, phi, axis, perm);
}

int main ()
{
  for (unsigned axis=0; axis < 3; axis++)
  {
    cerr << "Testing rotations about " << name[axis] << " axis" << endl;
    test_loop (axis, 0);
  }

  cerr << "All tests passed" << endl;
  return 0;
}

