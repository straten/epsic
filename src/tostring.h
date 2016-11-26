//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2004 - 2016 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

// psrchive/Util/units/tostring.h

#ifndef __TOSTRING_H
#define __TOSTRING_H

#include <string>
#include <sstream>
#include <limits>
#include <typeinfo>

// works around the circular dependence between Error and tostring
void raise (const char* name, const std::string& exception);

class ToString
{
  typedef std::ios_base::fmtflags fmtflags;

  unsigned precision;
  bool precision_set;

  fmtflags setf; 
  bool setf_set;

  fmtflags unsetf; 
  bool unsetf_set;

public:

  ToString () { reset_modifiers(); }

  void reset_modifiers () { precision_set = setf_set = unsetf_set = false; }

  void set_precision (unsigned p) { precision = p; precision_set = true; }
  void set_setf (fmtflags f) { setf = f; setf_set = true; }
  void set_unsetf (fmtflags f) { unsetf = f; unsetf_set = true; }

  template<class T>
  std::string operator () (const T& input) const
  {
    std::ostringstream ost;

    if (setf_set)
      ost.setf (setf);

    if (unsetf_set)
      ost.unsetf (unsetf);
    
    if (precision_set)
      ost.precision (precision);
    else
      ost.precision (std::numeric_limits<T>::digits10);
    
    ost << input;

    if (ost.fail())
      raise ("tostring", 
	     "failed to convert "+ std::string(typeid(T).name()) +" to string");

    return ost.str();
  }
};


/* the following global variables are not nested-call or multi-thread
   safe and should be used only when it is extremely difficult to pass
   the relevant arguments directly to the tostring function */

extern unsigned tostring_precision;
extern std::ios_base::fmtflags tostring_setf; 
extern std::ios_base::fmtflags tostring_unsetf; 

#define FMTFLAGS_ZERO std::ios_base::fmtflags(0)

template<class T>
std::string tostring (const T& input,
		      unsigned precision = std::numeric_limits<T>::digits10,
		      std::ios_base::fmtflags set = FMTFLAGS_ZERO,
		      std::ios_base::fmtflags unset = FMTFLAGS_ZERO)
{
  ToString tostr;

  if (tostring_setf)
    tostr.set_setf (tostring_setf);
  else if (set)
    tostr.set_setf (set);
  
  if (tostring_unsetf)
    tostr.set_unsetf (tostring_unsetf);
  else if (unset)
    tostr.set_unsetf (unset);
  
  if (tostring_precision)
    tostr.set_precision (tostring_precision);
  else
    tostr.set_precision (precision);

  return tostr( input );
}

template<class T>
T fromstring (const std::string& input)
{
  std::istringstream ist;

  ist.clear();
  ist.str(input);

  T retval;
  ist >> retval;

  if (ist.fail())
    raise ("fromstring", "failed to parse '"+ input +"'");

  return retval;
}

// string class specialization
inline std::string tostring (const std::string& input)
{
  return input;
}

// string class specialization
template<>
inline std::string fromstring<std::string> (const std::string& input)
{
  return input;
}

// char* specialization
inline std::string tostring (const char* input)
{
  return input;
}

#include "Error.h"

/*
  If you've already written a function that converts string to Type,
  and this function throws an exception of type Error when the string
  cannot be parsed, then you can use this template to implement an
  extraction operator (operator >>).  See Util/genutil/Types.C for an
  example. 
*/
template<typename Type, typename String2Type>
std::istream& extraction (std::istream& is, Type& t, String2Type string2type)
{
  std::streampos pos = is.tellg();
  std::string ss;

  try {
    is >> ss;
    t = string2type (ss);
  }
  catch (Error& e) {
    is.setstate(std::istream::failbit);
    is.seekg(pos);
  }
  return is;
}

#endif
