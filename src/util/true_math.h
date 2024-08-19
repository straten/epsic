/***************************************************************************
 *
 *   Copyright (C) 2024 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#ifndef __epsic_src_util_signbit_h
#define __epsic_src_util_signbit_h

#ifdef __cplusplus
extern "C" {
#endif

/*! works around ffast-math disabling signbit */
int true_signbit_float (float x);
int true_signbit_double (double x);
int true_signbit_long_double (long double x);

/*! works around fast-math disabling isfinite */
int true_finite_float (float x);
int true_finite_double (double x);
int true_finite_long_double (long double x);

#ifdef __cplusplus

} // extern "C"

// namespace for C++ wrappers of C functions
namespace true_math
{
  //! Overloaded wrappers of C functions
  inline int signbit (float x) { return true_signbit_float(x); }
  inline int signbit (double x) { return true_signbit_double(x); }
  inline int signbit (long double x) { return true_signbit_long_double(x); }

  //! Overloaded wrappers of C functions
  inline int finite (float x) { return true_finite_float(x); }
  inline int finite (double x) { return true_finite_double(x); }
  inline int finite (long double x) { return true_finite_long_double(x); }
}

#endif

#endif // __epsic_src_util_signbit_h

