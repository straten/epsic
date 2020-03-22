/***************************************************************************
 *
 *   Copyright (C) 2020 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "covariant.h"

double covariant_mode::modulation ()
{
  if (amps.size() == 0)
    coordinator->get();

  double retval = amps.front();
  amps.pop();
  return retval;
}

covariant_coordinator::covariant_coordinator (double covariance)
{
  a_in = b_in = 0;
  a_out = b_out = 0;

  covar[0][0] = covar[1][1] = 1.0;
  covar[0][1] = covar[1][0] = covariance;
}

void covariant_coordinator::set_modeA_input (modulated_mode* m)
{
  a_in = m;
}

void covariant_coordinator::set_modeB_input (modulated_mode* m)
{
  b_in = m;
}

modulated_mode* covariant_coordinator::get_modeA_output ()
{
  if (!a_out)
  {
    if (!a_in)
      throw std::runtime_error( "covariant_coordinator::get_modeA_output "
                                "modeA_input not set" );

    a_out = new covariant_mode ( a_in->get_source() );
    a_out->coordinator = this;
  }

  return a_out;
}

modulated_mode* covariant_coordinator::get_modeB_output ()
{
  if (!b_out)
  {
    if (!b_in)
      throw std::runtime_error( "covariant_coordinator::get_modeB_output "
                                "modeB_input not set" );

    b_out = new covariant_mode ( b_in->get_source() );
    b_out->coordinator = this;
  }

  return b_out;
}

void covariant_coordinator::get()
{
  if (!a_in)
    throw std::runtime_error( "covariant_coordinator::get "
                              "modeA_input not set" );
  if (!b_in)
    throw std::runtime_error( "covariant_coordinator::get "
                              "modeB_input not set" );
  if (!a_out)
    throw std::runtime_error( "covariant_coordinator::get "
                              "modeA_output not set" );
  if (!b_out)
    throw std::runtime_error( "covariant_coordinator::get "
                              "modeB_output not set" );

  Vector<2,double> amps (a_in->modulation(), b_in->modulation());
  Vector<2,double> result = covar * amps;

  a_out->amps.push( result[0] );
  b_out->amps.push( result[1] );
}

