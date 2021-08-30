/***************************************************************************
 *
 *   Copyright (C) 2021 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "modulated.h"

using namespace std;

square_modulated_mode::square_modulated_mode (modulated_mode* _mode,
					      unsigned _width,
					      unsigned sample_size)
  : modulated_mode( _mode->get_source() )
{
  width = _width;
  current = _width;
  mod = _mode;

  compute_cross_correlation (sample_size);
}

void square_modulated_mode::compute_cross_correlation (unsigned sample_size)
{
  cross_correlation.resize( width );

  if (sample_size <= width)
  {
    for (unsigned i=0; i<width; i++)
      cross_correlation[i] = 1.0;

    return;
  }

  vector<double> matrix ( sample_size * sample_size, 0.0 );

  unsigned cur = 0;
  unsigned populations = 0;
  while (cur != width)
  {
    unsigned ioff=0;
    while (ioff < sample_size)
    {
      unsigned end = width;
      if (ioff + end > sample_size)
	end = sample_size - ioff;

      if (cur == width)
	cur = 0;

#if _DEBUG
      cerr << "ioff=" << ioff << " cur=" << cur << " end=" << end << endl;
#endif
      
      for (; cur < end ; cur++)
      {
	for (unsigned col=cur; col < end; col++)
	  matrix[ ioff*sample_size + (ioff+col-cur) ] += 1.0;
	ioff ++;
      }
    }

#if _DEBUG
    for (unsigned irow=0; irow < sample_size; irow++)
    {
      for (unsigned icol=0; icol < sample_size; icol++)
	{
	  unsigned jrow=irow;
	  unsigned jcol=icol;
	  if (jcol < jrow)
	    swap (jcol, jrow);
	  cerr << matrix[ jrow*sample_size + jcol ] << " ";
	}
      cerr << endl;
    }
#endif
    
    populations ++;
  }

#if _DEBUG
  cerr << "square_modulated_mode::compute_cross_correlation "
       << populations << " populations" << endl;
#endif
  
  for (unsigned ilag=0; ilag < width; ilag++)
  {
    unsigned nrow = sample_size - ilag;

    double sum = 0;
    for (unsigned irow=0; irow < nrow; irow++)
      sum += matrix[ irow*sample_size + (irow+ilag) ];
    
    cross_correlation[ilag] = sum / (nrow * populations);
  }
}
