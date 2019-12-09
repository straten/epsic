/***************************************************************************
 *
 *   Copyright (C) 2016 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "mode.h"
#include "Quaternion.h"
#include "Pauli.h"

using namespace std;

mode::mode ()
{
  normal = 0;
  rms = 0.5;

  set_Stokes (Stokes<double>(1.0));
}

void mode::set_Stokes (const Stokes<double>& _mean)
{
  mean = _mean;

  Quaternion<double,Hermitian> root = sqrt (natural(mean));
  polarizer = convert (root);
}

Spinor<double> mode::get_field ()
{
  if (!normal)
    throw std::runtime_error( "mode::get_field - BoxMuller not set");

  BoxMuller& gasdev = *normal;

  complex<double> x (rms * gasdev(), rms * gasdev());
  complex<double> y (rms * gasdev(), rms * gasdev());

  Spinor<double> e (x, y);

  return polarizer * e;
}

#if 0
void mode::powerlaw (Spinor<double>& e)
{
  double y = random_double(); 
  double index = -4;

  e *= pow( pow(inner_scale,index+1) * (1.0 - y), 1.0/(index+1) );
}
#endif

