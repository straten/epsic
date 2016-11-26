/***************************************************************************
 *
 *   Copyright (C) 2004 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
#include "Error.h"

#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

using namespace std;

bool Error::verbose = false;
bool Error::complete_abort = false;

Error::Error ()
{
  code = Undefined;
}

Error::Error (ErrorCode c, string func, const char* msg, ...)
{
  static char buf[1024];
  char* this_msg = 0;

  if (msg) {
    va_list args;
  
    va_start(args, msg);
    vsnprintf(buf, 1024, msg, args);
    va_end(args);

    this_msg = buf;
  }
  construct (c, func, this_msg);
}


//! Error with string message
Error::Error (ErrorCode c, string func, string msg)
{
  construct (c, func, msg.c_str());
}

void Error::construct (ErrorCode c, const string& func, const char* msg)
{
  code = c; 

  if (!func.empty())
    functions.push_back(func);

  if (msg)
    message = msg;

  errno_check();

  if (verbose)
    cerr << "Error::construct" << *this << endl;

  if (complete_abort)
    abort ();
}

void Error::errno_check ()
{
  if (code == FailedSys && errno != 0) {
    message += " - ";
    message += strerror (errno);
  }
}

Error::~Error ()
{
}

//! Add function name to the list
const Error& Error::operator += (const char* func)
{
  functions.push_back (func);
  if (verbose)
    report (cerr);
  return *this;
}
  
//! Add function name to the list
const Error& Error::operator += (const string& func)
{
  functions.push_back (func);
  if (verbose)
    report (cerr);
  return *this;
}

void Error::report (ostream& ostr) const
{
  ostr << endl;

  // print the list of functions in reverse order
  if (!functions.empty()) {
    ostr << "Error::stack" << endl;
    for (int i = int(functions.size())-1; i>=0; i--)
      ostr << "\t" << functions[i] << endl;
  }

  ostr << "Error::" << err2str (code) << endl;

  // print the error message, if any
  if (!message.empty())
    ostr << "Error::message\n\t" << message << endl;
}

string Error::warning () const
{
  string function;

  if (!functions.empty())
    function = functions.back();

  return function + " WARNING: " + message;
}


const string Error::get_message() const
{
  return message;
}

ostream& operator<< (ostream& ostr, const Error& error)
{
  error.report (ostr);
  return ostr;
}

const char* Error::err2str (ErrorCode code)
{
#define CASESTR(x) case x: return #x
  switch (code)  {
    CASESTR(Undefined);
    CASESTR(BadAllocation);
    CASESTR(InvalidPointer);
    CASESTR(InvalidParam);
    CASESTR(InvalidState);
    CASESTR(InvalidRange);
    CASESTR(FileNotFound);
    CASESTR(EndOfFile);
    CASESTR(FailedCall);
    CASESTR(FailedSys);
  default: return "invalid";
  }
#undef CASESTR
}
