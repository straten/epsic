/***************************************************************************
 *
 *   Copyright (C) 2006 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "Convention.h"

#include <stdlib.h>

using namespace std;

ostream& operator << (ostream& os, Convention::Basis basis)
{
  switch (basis)
  {
  case Convention::Linear:
    return os << "lin";
  case Convention::Circular:
    return os << "cir";
  case Convention::Elliptical:
    return os << "ell";
  }
  return os;
}

istream& operator >> (istream& is, Convention::Basis& basis)
{
  std::streampos pos = is.tellg();
  string unit;
  is >> unit;

  if (unit == "lin" || unit == "Linear")
    basis = Convention::Linear;
  else if (unit == "cir" || unit == "circ" || unit == "Circular")
    basis = Convention::Circular;
  else if (unit == "ell" || unit == "Elliptical")
    basis = Convention::Elliptical;
  else {

    // replace the text and try to parse a number
    is.seekg (pos);

    int code = -1;
    is >> code;

    if (!is.fail())
    {
      switch ((Convention::Basis)code)
      {
      case Convention::Linear:
	basis = Convention::Linear; break;
      case Convention::Circular:
	basis = Convention::Circular; break;
      case Convention::Elliptical:
	basis = Convention::Elliptical; break;
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

ostream& operator << (ostream& os, Convention::Hand hand)
{
  return output (os, hand);
}

istream& operator >> (istream& is, Convention::Hand& hand)
{
  return input (is, hand);
}

ostream& operator << (ostream& os, Convention::Argument arg)
{
  return output (os, arg);
}

istream& operator >> (istream& is, Convention::Argument& arg)
{
  return input (is, arg);
}

