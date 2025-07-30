//-*-C++-*-

/***************************************************************************
 *
 *   Copyright (C) 2025 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// thrust code needs to be compiled by nvcc ...
extern int test_cuSpinor();

// ... but we don't have rules in place to link executables using nvcc
int main()
{
  return test_cuSpinor();
  return 0;
}
