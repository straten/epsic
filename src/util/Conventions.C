/***************************************************************************
 *
 *   Copyright (C) 2006 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "Conventions.h"

#include <stdlib.h>

using namespace std;

ostream& operator << (ostream& os, Signal::Basis basis)
{
  switch (basis)
  {
  case Signal::Linear:
    return os << "lin";
  case Signal::Circular:
    return os << "cir";
  case Signal::Elliptical:
    return os << "ell";
  }
  return os;
}

istream& operator >> (istream& is, Signal::Basis& basis)
{
  std::streampos pos = is.tellg();
  string unit;
  is >> unit;

  if (unit == "lin" || unit == "Linear")
    basis = Signal::Linear;
  else if (unit == "cir" || unit == "circ" || unit == "Circular")
    basis = Signal::Circular;
  else if (unit == "ell" || unit == "Elliptical")
    basis = Signal::Elliptical;
  else {

    // replace the text and try to parse a number
    is.seekg (pos);

    int code = -1;
    is >> code;

    if (!is.fail())
    {
      switch ((Signal::Basis)code)
      {
      case Signal::Linear:
	basis = Signal::Linear; break;
      case Signal::Circular:
	basis = Signal::Circular; break;
      case Signal::Elliptical:
	basis = Signal::Elliptical; break;
      }
    }

  }

  return is;
}

template<typename T>
ostream& output (ostream& os, T argument)
{
  ostream::fmtflags flags = os.setf (ostream::showpos);
  os << (int) argument;
  os.flags (flags);
  return os;
}

template<typename T>
istream& input (istream& is, T& argument)
{
  int code = 0;
  is >> code;

  // Hand and Argument may equal only +/- 1
  if (abs(code) != 1)
    is.setstate (ios::failbit);

  argument = (T) code;
  return is;
}

ostream& operator << (ostream& os, Signal::Hand hand)
{
  return output (os, hand);
}

istream& operator >> (istream& is, Signal::Hand& hand)
{
  return input (is, hand);
}

ostream& operator << (ostream& os, Signal::Argument arg)
{
  return output (os, arg);
}

istream& operator >> (istream& is, Signal::Argument& arg)
{
  return input (is, arg);
}

