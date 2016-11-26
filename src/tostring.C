/***************************************************************************
 *
 *   Copyright (C) 2006-2010 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "tostring.h"

unsigned tostring_precision (0);
std::ios_base::fmtflags tostring_setf = FMTFLAGS_ZERO;
std::ios_base::fmtflags tostring_unsetf = FMTFLAGS_ZERO;

/*
  operator << (Error&, T) uses tostring, and tostring throws an Error

  This was never a problem until gcc 4.6, which will no longer
  postpone template interpretation until all of the required symbols
  (classes, functions, etc.) have been declared.

  Therefore, this silly wrapper has to be written to explicitly
  postpone dependency on the Error class.
*/
void raise (const char* name, const std::string& exception)
{
  throw Error (InvalidState, name, exception);
}
