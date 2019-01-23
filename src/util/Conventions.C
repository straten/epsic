/***************************************************************************
 *
 *   Copyright (C) 2006 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "Conventions.h"

#include <stdlib.h>

using namespace std;

ostream& Signal::operator << (ostream& os, Basis basis)
{
  switch (basis) {
  case Linear:
    return os << "lin";
  case Circular:
    return os << "cir";
  case Elliptical:
    return os << "ell";
  }
  return os;
}

istream& Signal::operator >> (istream& is, Basis& basis)
{
  std::streampos pos = is.tellg();
  string unit;
  is >> unit;

  if (unit == "lin" || unit == "Linear")
    basis = Linear;
  else if (unit == "cir" || unit == "circ" || unit == "Circular")
    basis = Circular;
  else if (unit == "ell" || unit == "Elliptical")
    basis = Elliptical;
  else {

    // replace the text and try to parse a number
    is.seekg (pos);

    int code = -1;
    is >> code;

    if (!is.fail()) {
      switch ((Basis)code) {
      case Linear:
	basis = Linear; break;
      case Circular:
	basis = Circular; break;
      case Elliptical:
	basis = Elliptical; break;
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

ostream& Signal::operator << (ostream& os, Hand hand)
{
  return output (os, hand);
}

istream& Signal::operator >> (istream& is, Hand& hand)
{
  return input (is, hand);
}

ostream& Signal::operator << (ostream& os, Argument arg)
{
  return output (os, arg);
}

istream& Signal::operator >> (istream& is, Argument& arg)
{
  return input (is, arg);
}
