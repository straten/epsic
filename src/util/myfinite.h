/***************************************************************************
 *
 *   Copyright (C) 2023 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#ifndef __myfinite_H
#define __myfinite_H

#ifdef __cplusplus
extern "C" {
#endif

/*! works around fast-math disabling isfinite */
int myfinite (double x);

#ifdef __cplusplus
}
#endif

#endif

