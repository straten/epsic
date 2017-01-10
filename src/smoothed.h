//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2016 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// epsic/src/util/smoothed.h

#ifndef __smoothed_H
#define __smoothed_H

#include "mode.h"

/***************************************************************************
 *
 *  a boxcar-smoothed mode of electromagnetic radiation
 *
 ***************************************************************************/

class boxcar_mode : public mode_decorator
{
  std::vector< Spinor<double> > instances;
  unsigned smooth;
  unsigned current;

  void setup()
  {
    current = 0;
    instances.resize (smooth);
    for (unsigned i=1; i<smooth; i++)
      instances[i] = source->get_field();
  }

public:

  boxcar_mode (mode* s, unsigned n) : mode_decorator(s) { smooth = n; }

  Spinor<double> get_field ()
  {
    if (instances.size() < smooth)
      setup ();

    instances[current] = source->get_field();
    current = (current + 1) % smooth;

    Spinor<double> result;
    for (unsigned i=0; i<smooth; i++)
      result += instances[i];

    result /= sqrt(smooth);
    return result;
  }
};

#endif
