//-*-C++-*-
/***************************************************************************
 *
 *   Copyright (C) 2004 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
/* $Source: /cvsroot/psrchive/psrchive/Util/units/complex_promote.h,v $
   $Revision: 1.6 $
   $Date: 2006/10/06 21:13:55 $
   $Author: straten $ */

#ifndef __complex_promote_H
#define __complex_promote_H

#include <complex>
#include "PromoteTraits.h"

#ifdef PROMOTE_TRAITS_SPECIALIZE

template <class T, class U>
class PromoteTraits< std::complex<T>, std::complex<U> >
{
  public:
    typedef std::complex<typename PromoteTraits<T,U>::promote_type> 
            promote_type;
};

#endif

//! std::complex addition
template<typename T, typename U>
const std::complex<typename PromoteTraits<T,U>::promote_type>
operator + (const std::complex<T>& a, const std::complex<U>& b)
{
  std::complex<typename PromoteTraits<T,U>::promote_type> ret(a);
  ret+=b;
  return ret;
}

//! std::complex subtraction
template<typename T, typename U>
const std::complex<typename PromoteTraits<T,U>::promote_type>
operator - (const std::complex<T>& a, const std::complex<U>& b)
{
  std::complex<typename PromoteTraits<T,U>::promote_type> ret(a);
  ret-=b;
  return ret;
}

//! std::complex multiplication
template<typename T, typename U>
const std::complex<typename PromoteTraits<T,U>::promote_type>
operator * (const std::complex<T>& a, const std::complex<U>& b)
{
  std::complex<typename PromoteTraits<T,U>::promote_type> ret(a);
  ret*=b;
  return ret;
}

#endif

