/***************************************************************************
 *
 *   Copyright (C) 2016 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "mode.h"
#include "Quaternion.h"
#include "Pauli.h"

epsic::mode::mode ()
{
  normal = 0;
  rms = 0.5;

  set_Stokes (Stokes<double>(1.0));
}

void epsic::mode::set_Stokes (const Stokes<double>& _mean)
{
  mean = _mean;

  Quaternion<double,Hermitian> root = sqrt (natural(mean));
  polarizer = convert (root);
}

Spinor<double> epsic::mode::get_field ()
{
  if (!normal)
    throw std::runtime_error( "epsic::mode::get_field - BoxMuller not set");

  BoxMuller& gasdev = *normal;

  std::complex<double> x (rms * gasdev(), rms * gasdev());
  std::complex<double> y (rms * gasdev(), rms * gasdev());

  Spinor<double> e (x, y);

  return polarizer * e;
}


