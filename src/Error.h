//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2004 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

/* $Source: /cvsroot/psrchive/psrchive/Util/units/Error.h,v $
   $Revision: 1.9 $
   $Date: 2008/01/02 22:38:05 $
   $Author: straten $ */

#ifndef __Error_h
#define __Error_h

#include <iostream>
#include <vector>
#include <string>

//! List of error types
enum ErrorCode {
  //! undefined error
  Undefined,
  //! operator new has returned null
  BadAllocation,
  //! General bad pointer
  InvalidPointer,
  //! general bad parameter
  InvalidParam,
  //! general invalid state
  InvalidState,
  //! index out of range, or range mismatch
  InvalidRange,
  //! file not found
  FileNotFound,
  //! end of file
  EndOfFile,
  //! failure of a function called from another lib
  FailedCall,
  //! failure of a system call that sets errno
  FailedSys
};

//! A convenient exception handling class
/*
  It is useful to throw an exception with both an error code and a
  text message.  As well, on catching an exception, it is useful if a
  function or method can add its name to an "exception stack",
  enabling the error to be traced back to its origin.  The static
  method, Error::handle_signals(), may also be invoked so that signal
  interrupts, such as SIGFPE and SIGSEGV, will be intercepted and an
  Error exception will be thrown. 
*/
class Error {

  public:

  //! convert an error code to a character string
  static const char* err2str (ErrorCode code);

  //! During constructor, print messages
  static bool verbose;

  //! After construction raise a sigkill
  static bool complete_abort;

  //! Error with optional printf-style message
  Error (ErrorCode c, std::string func, const char* msg=0, ...);

  //! Error with string message
  Error (ErrorCode c, std::string func, std::string msg);
  
  //! Virtual destructor enables inheritance
  virtual ~Error ();

  //! Add function name to the list
  const Error& operator += (const char* func);
  
  //! Add function name to the list
  const Error& operator += (const std::string& func);
  
  //! Write error and auxilliary info to ostream
  virtual void report (std::ostream& ostr) const;
  
  //! Write function name and error message to screen
  virtual std::string warning () const;

  //! Get the error message
  const std::string get_message() const;

  //! Get the error code
  ErrorCode get_code () const { return code; }

  //! Add to the Error message
  template<class T>
  friend Error& operator<< (Error& error, const T& t);

protected:

  //! the error code
  ErrorCode code;
  //! auxilliary message
  std::string message;
  //! function name stack
  std::vector<std::string> functions;

  //! Null constructor is protected
  Error ();

  //! If errno is set, adds a string to message
  void errno_check ();
  
  //! Contruct the basic deal
  void construct (ErrorCode c, const std::string& func, const char* msg);
};

//! Convenience output operator in global namespace
std::ostream& operator<< (std::ostream& ostr, const Error& error);

#include "tostring.h"

template<class T>
Error& operator<< (Error& error, const T& t)
{
  error.message += tostring(t);
  return error;
}

#endif
