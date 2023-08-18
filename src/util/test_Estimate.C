/***************************************************************************
 *
 *   Copyright (C) 2003 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
#include "Estimate.h"
#include "Stokes.h"

#include <complex>

using namespace std;

int main ()
{
  // ///////////////////////////////////////////////////////////////////////
  //
  // test addition
  //
  
  Estimate<double> e1 (5, 1);
  Estimate<double> e2 (7, 1);
  
  Estimate<double> result = e1 + e2;
  
  if (result.val != 12.0 || result.var != 2.0) {
    cerr << "test_Estimate: Estimate addition error " << result << endl;
    return -1;
  }
      
  
  Estimate<double> x (4.6, 0.09);
  Estimate<double> y (-0.5, 0.01);

  // ///////////////////////////////////////////////////////////////////////
  //
  // test multiplication
  //

  result = x*y;

  // dr/dx = y = -0.5
  // dr/dy = x = 4.6
  // r.var = 0.25 * 0.09 + 21.16 * 0.01 = 0.2341

  if (fabs(result.val - -2.3) > 1e-6 || fabs(result.var - 0.2341) > 1e-6) {
    cerr << "test_Estimate: Estimate multiplication error " << result << endl;
    return -1;
  }

  // ///////////////////////////////////////////////////////////////////////
  //
  // test division
  //

  result = x/y;

  // dr/dx = 1/y = -2
  // dr/dy = -x/y^2 = -18.4
  // sigma_r^2 = 4 * 0.09 + 338.56 * 0.01 = 3.7456

  if (fabs(result.val - -9.2) > 1e-6 || fabs(result.var - 3.7456) > 1e-6) {
    cerr << "test_Estimate: Estimate division error " << result << endl;
    return -1;
  }

  // ///////////////////////////////////////////////////////////////////////
  //
  // test the accumulation of mean 
  //

  MeanEstimate<double> mean;

  mean += e1;
  mean += e2;

  result = mean.get_Estimate();

  cerr << "test_Estimate: mean of e1=" << e1 << " and e2=" << e2 <<
    "\n  = " << result << endl;

  if (result.val != 6.0 || result.var != 0.5) {
    cerr << "test_Estimate: MeanEstimate error " << result << endl;
    return -1;
  }

  // ///////////////////////////////////////////////////////////////////////
  //
  // test the accumulation of mean in radians
  //

  Estimate<double> rad1 (M_PI/4.0, 0.01);
  Estimate<double> rad2 (2.0*M_PI-M_PI/4.0, 0.01);

  MeanRadian<double> radmean;

  radmean += rad1;
  radmean += rad2;

  Estimate<double> rad3 (radmean);

  cerr << "test_Estimate: mean of rad1=" << rad1 << " and rad2=" << rad2 <<
    "\n  = " << rad3 << endl;
    
  if (fabs(rad3.val) > 1e-10 || fabs(rad3.var - 0.005) > 1e-10) {
    cerr << "test_Estimate: MeanRadian error " << rad3 << endl;
    return -1;
  }

  Estimate<double> StokesI (0.967525,0.000659444);
  Estimate<double> StokesQ (0.202817,0.000658172);
  Estimate<double> StokesU (0.0604878,0.000660645);
  Estimate<double> StokesV (-0.0498089,0.00068193);

  Estimate<double> inv (0.888829,0.00129698);

  cerr << "inv=" << inv << endl;
  inv = sqrt(inv);
  cerr << "S=sqrt(inv)=" << inv << endl;

  cerr << "I/S=" << StokesI/inv << endl;
  cerr << "Q/S=" << StokesQ/inv << endl;


  // norm=[h:(1.02625±0.00117418),(0.215127±0.000759488),(0.0641592±0.000744965),(-0.0528321±0.000768368)]

  // ///////////////////////////////////////////////////////////////////////
  //
  // test complex<Estimate> operations
  //

  std::complex< Estimate<double> > a (Estimate<double> (1.2,0.1), Estimate<double> (3.4, 0.1));
  std::complex< Estimate<double> > b (Estimate<double> (5.6,0.1), Estimate<double> (-7.8, 0.1));

  std::complex< Estimate<double> > cresult;

  const double tol = 1e-10;

  cerr << "Complex values:" << endl;

  cresult = a*b;
  cerr << a << " * " << b << " = " << cresult << endl;
  if (fabs(cresult.real().val-33.24)>tol || fabs(cresult.imag().val-9.68)>tol)
  {
    cerr << "test_Estimate: complex multiply error " << cresult << endl;
    return -1;
  }

  result = norm(a);
  cerr << "norm(" << a << ") = " << result << endl;
  if (fabs(result.val-13.0)>tol)
  {
    cerr << "test_Estimate: complex norm error " << cresult << endl;
    return -1;
  }

  
  return 0;
}

