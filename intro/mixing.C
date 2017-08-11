
#include "Plot3D.h"
#include "MEAL/Gaussian.h"

#include <cpgplot.h>

#include <iostream>
using namespace std;

#include <math.h>
#include <unistd.h>

void drawAxis (pgplot::Plot3D& volume, const Cartesian& origin,
	       double xlen, double ylen, double zlen, bool arrows=false)
{
  Cartesian x0 (1,0,0);
  Cartesian y0 (0,1,0);
  Cartesian z0 (0,0,1);

  if (arrows) {
    volume.arrow (origin-xlen*x0, origin+xlen*x0);
    volume.arrow (origin-ylen*y0, origin+ylen*y0);
    volume.arrow (origin-zlen*z0, origin+zlen*z0);
  }
  else {
    volume.move (origin-xlen*x0); volume.draw ( origin+xlen*x0);
    volume.move (origin-ylen*y0); volume.draw ( origin+ylen*y0);
    volume.move (origin-zlen*z0); volume.draw ( origin+zlen*z0);
  }

}

void drawShape (pgplot::Plot3D& volume, const Cartesian& origin,
		const Cartesian& ordaxis, const Cartesian& absaxis,
		vector<double>& shape, double ord1, double ord2, int fill=0,
		unsigned ordi1=0, unsigned ordi2=0)
{
  if (ordi2==0)
    ordi2 = shape.size();

  double step = 0;

  if (shape.size() > 1)
    step = (ord2 - ord1) / (shape.size()-1);
  
  vector<Cartesian> poly;

  // help to make the shape when function does not meet axis
  Cartesian start = origin + (ord1 + step*double(ordi1)) * ordaxis;
  Cartesian end = origin + (ord1 + step*double(ordi2-1)) * ordaxis;


  if (fill) {
    int fs=0;
    cpgqfs(&fs);

    if (fs == 3)
      // make hatching run parallel to abscissa axis
      volume.set_hatch (absaxis);

    // this zeros the edges
    poly.push_back(start);
  }

  volume.move(start);

  for (unsigned iord=ordi1; iord<ordi2; iord++) {
    Cartesian pt =
      origin
      + (ord1 + step*double(iord)) * ordaxis
      + shape[iord] * absaxis;

    if (fill)
      poly.push_back(pt);

    if (!fill || fill==2)
      volume.draw(pt);

  }

  if (!fill || fill==2)
    volume.draw(end);

  if (fill) {
    // this zeros the edges
    poly.push_back(end);
    volume.poly (poly);
  }

}


int main (int argc, char** argv)
{
  // position angles defining point of view
  float theta = 30; 
  float phi = 30;

  // fractional width of passband kept in band-limited signal
  float bandwidth = 0.2;

  int gotc;

  bool quadrature = true;

  string device = "/xs";

  while ((gotc = getopt(argc, argv, "hb:D:p:ult:v")) != -1) {
    switch (gotc) {

    case 'h':
      cerr << endl <<
	" -b bandwidth\n"
	" -D device\n"
	" -t theta view\n"
	" -p phi view\n"
	" -u upper single-sideband\n"
	" -l lower single-sideband\n"
	   << endl;
      return 0;

    case 'b':
      bandwidth = atof (optarg);
      break;

    case 'D':
      device = optarg;
      break;

    case 't':
      theta = atof (optarg);
      break;
    case 'p':
      phi = atof (optarg);
      break;

    case 'u':
      quadrature = false;
      break;

    case 'v':
      pgplot::Plot3D::verbose = true;
      break;
    }
  }

  float width = 2.5;
  float height = 2.1;

  cpgbeg  (0, device.c_str(), 0, 0);
  cpgpap (0.0, height/width);
  cpgswin (-width,width,-height,height);
  cpgslw (3);

  // ////////////////////////////////////////////////////////////////
  //
  // make the groovy shapes to be plotted
  //
  // ////////////////////////////////////////////////////////////////

  MEAL::Gaussian hump1;
  hump1.set_centre (.65);
  hump1.set_height (.30);
  hump1.set_width (.14);
  
  MEAL::Gaussian hump2;
  hump2.set_centre (.37);
  hump2.set_height (.12);
  hump2.set_width (.11);
  
  int npts = 50;
  vector<double> Xw (npts);
  vector<double> Sum_w (npts);
  vector<double> Dif_w (npts);

  for (int iw=0; iw<npts; iw++)
  {
    double x = double(iw)/npts;
    hump1.set_abscissa ( x );
    hump2.set_abscissa ( x );
    
    Xw[iw] = hump1.evaluate() + hump2.evaluate();
  }
  
  for (int iw=0; iw<npts; iw++)
  {
    Sum_w[iw] = Xw[iw] + Xw[npts-iw-1];
    Dif_w[iw] = Xw[iw] - Xw[npts-iw-1];
  }

  // ////////////////////////////////////////////////////////////////
  //
  // set up the origins and axis
  //
  // ////////////////////////////////////////////////////////////////

  pgplot::Plot3D volume;
  //volume.set_frame (Cartesian(0.0,0.0), Cartesian(1,1));
  volume.set_camera (theta,phi);

  float dist = 1.0;
  Cartesian hsep = dist*volume.get_yaxis ();
  Cartesian wsep = volume.get_xaxis ();

  // original spectrum
  Cartesian Ox = 1.5 * hsep;

  // //////////////////////////////////////////////////
  // these apply to quadrature mode only

  // in-phase spectrum
  Cartesian Oi = .5 * hsep - wsep;
  // quadrature spectrum
  Cartesian Oq = .5 * hsep + wsep;

  // band-limited, in-phase spectrum
  Cartesian Obi = -.5 * hsep - wsep;
  // band-limited, quadrature spectrum
  Cartesian Obq = -.5 * hsep + wsep;

  // //////////////////////////////////////////////////
  // these apply to single-sideband mode only

  // band-limited, original spectrum
  Cartesian Ob = .5 * hsep;

  // upper-sideband signal
  Cartesian Ou = -.5 * hsep;

  // band-limited, baseband spectrum
  Cartesian Obx = -1.5 * hsep;

  // the three dimensions: frequency (nu), real (Re), and imaginary (Im)
  Cartesian nu (0,1,0);
  Cartesian Re (0,0,1);
  Cartesian Im (1,0,0);

  // set the colours to help differentiate
  int cRe = 16;
  int cIm = 17;
  float fRe = 0.45;
  float fIm = 0.75;
  cpgscr (cRe, fRe, fRe, fRe);
  cpgscr (cIm, fIm, fIm, fIm);

  // fill style: hatching
  cpgsfs(3);

  // ////////////////////////////////////////////////////////////////
  //
  // top plot, X(w)
  //
  // ////////////////////////////////////////////////////////////////

  cpgsci(cIm);
  drawShape (volume, Ox, nu, Im, Xw, 0, 1, 2);
  drawShape (volume, Ox, nu, -Im, Xw, 0, -1, 2);
		     
  cpgsci(cRe);	     
  drawShape (volume, Ox, nu, Re, Xw, 0, 1, 2);
  drawShape (volume, Ox, nu, Re, Xw, 0, -1, 2);

  if (quadrature) {

    // ////////////////////////////////////////////////////////////////
    //
    // in phase, I(w)
    //
    // ////////////////////////////////////////////////////////////////
    cpgsci(cIm);
    drawShape (volume, Oi, nu, Im, Xw, .5, 1.5, 2);
    drawShape (volume, Oi, nu, -Im, Xw, -.5, -1.5, 2);
    drawShape (volume, Oi, nu, Im, Dif_w, -.5, .5, 2);
    
    // show the unaliased forms as dotted lines
    // cpgsls(2);
    // drawShape (volume, Oi, nu, -Im, Xw, .5, -.5);
    // drawShape (volume, Oi, nu, Im, Xw, -.5, .5);
    // cpgsls(1);
    
    cpgsci(cRe);
    drawShape (volume, Oi, nu, Re, Xw, .5, 1.5, 2);
    drawShape (volume, Oi, nu, Re, Xw, -.5, -1.5, 2);
    drawShape (volume, Oi, nu, Re, Sum_w, -.5, .5, 2);
    
    // cpgsls(2);
    // drawShape (volume, Oi, nu, Re, Xw, .5, -.5);
    // drawShape (volume, Oi, nu, Re, Xw, -.5, .5);
    // cpgsls(1);
    
    // ////////////////////////////////////////////////////////////////
    //
    // quadrature, Q(w) 
    //
    // ////////////////////////////////////////////////////////////////
    
    // this one is behind
    cpgsci(cRe);
    drawShape (volume, Oq, nu, Im, Xw, .5, 1.5, 2);
    
    cpgsci(cIm);
    drawShape (volume, Oq, nu, -Re, Xw, .5, 1.5, 2);
    drawShape (volume, Oq, nu, -Re, Xw, -.5, -1.5, 2);
    drawShape (volume, Oq, nu, Re, Sum_w, -.5, .5, 2);
    
    cpgsci(cRe);
    drawShape (volume, Oq, nu, -Im, Xw, -.5, -1.5, 2);
    drawShape (volume, Oq, nu, -Im, Dif_w, -.5, .5, 2);
    
  }

  // ////////////////////////////////////////////////////////////////
  //
  // band-limited
  //
  // ////////////////////////////////////////////////////////////////

  int bpts = int (float(npts)*bandwidth);
  int bs = (npts - bpts)/2;
  int be = (npts + bpts)/2;

  if (quadrature) {
    cpgsci(cIm);
    drawShape (volume, Obi, nu, Im, Dif_w, -.5, .5, 2, bs,be);
    cpgsci(cRe);
    drawShape (volume, Obi, nu, Re, Sum_w, -.5, .5, 2, bs,be);
    
    cpgsci(cIm);
    drawShape (volume, Obq, nu, Re, Sum_w, -.5, .5, 2, bs,be);
    cpgsci(cRe);
    drawShape (volume, Obq, nu, -Im, Dif_w, -.5, .5, 2, bs,be);
    
    cpgsci(cIm);
    drawShape (volume, Obx, nu, Im, Xw, -.5, .5, 2, bs,be);
    cpgsci(cRe);
    drawShape (volume, Obx, nu, Re, Xw, -.5, .5, 2, bs,be);
  }
  else {

    // draw the band-limited passband
    cpgsci(cIm);
    drawShape (volume, Ob, nu, Im, Xw, 0, 1, 2, bs,be);
    drawShape (volume, Ob, nu, -Im, Xw, 0, -1, 2, bs,be);
    
    cpgsci(cRe);	     
    drawShape (volume, Ob, nu, Re, Xw, 0, 1, 2, bs,be);
    drawShape (volume, Ob, nu, Re, Xw, 0, -1, 2, bs,be);

    float lo = float(bs)/float(npts);

    // draw the band-limited upper-sideband passband
    cpgsci(cIm);
    drawShape (volume, Ou, nu, Im, Xw, -lo, 1-lo, 2, bs,be);
    drawShape (volume, Ou, nu, Im, Xw, lo, 1+lo, 2, bs,be);

    drawShape (volume, Ou, nu, -Im, Xw, lo, -1+lo, 2, bs,be);
     drawShape (volume, Ou, nu, -Im, Xw, -lo, -1-lo, 2, bs,be);
   
    cpgsci(cRe);	     
    drawShape (volume, Ou, nu, Re, Xw, -lo, 1-lo, 2, bs,be);
    drawShape (volume, Ou, nu, Re, Xw, lo, 1+lo, 2, bs,be);

    drawShape (volume, Ou, nu, Re, Xw, lo, -1+lo, 2, bs,be);
    drawShape (volume, Ou, nu, Re, Xw, -lo, -1-lo, 2, bs,be);

    // draw the band-limited baseband passband
    cpgsci(cIm);
    drawShape (volume, Obx, nu, Im, Xw, -lo, 1-lo, 2, bs,be);
    drawShape (volume, Obx, nu, -Im, Xw, lo, -1+lo, 2, bs,be);
    
    cpgsci(cRe);	     
    drawShape (volume, Obx, nu, Re, Xw, -lo, 1-lo, 2, bs,be);
    drawShape (volume, Obx, nu, Re, Xw, lo, -1+lo, 2, bs,be);

  }

  float labh = .45;

  Cartesian ql = .25*Re - .6 * bandwidth * nu;

  // draw all of the axis for the above plots
  cpgsci(1);
  drawAxis (volume, Ox, .25, 1.1, .25);
  volume.text (Ox + labh*Re, "X(\\gn)");

  if (quadrature) {

    drawAxis (volume, Oi, .25, 1.6, .25);
    volume.text (Oi + labh*Re, "I(\\gn)");

    drawAxis (volume, Oq, .25, 1.6, .25);
    volume.text (Oq + labh*Re, "Q(\\gn)");

    drawAxis (volume, Obi, .25, .5, .25);
    volume.text (Obi + ql, "I\\db\\u(\\gn)", 1);

    drawAxis (volume, Obq, .25, .5, .25);
    volume.text (Obq + ql, "Q\\db\\u(\\gn)", 1);

  }
  else {

    drawAxis (volume, Ob, .25, 1.1, .25);
    volume.text (Ob + labh*Re, "X\\db\\u(\\gn)");

    drawAxis (volume, Ou, .25, 1.6, .25);
    volume.text (Ou + labh*Re, "X\\dm\\u(\\gn)");
  }
  
  drawAxis (volume, Obx, .25, .5, .25);

  if (quadrature)
    volume.text (Obx + labh*Re, "Z\\db\\u(\\gn)");
  else
    volume.text (Obx + labh*Re, "X\\du\\u(\\gn)");

  cerr << "finished plotting" << endl;
  cpgend();
  cerr << "exiting" << endl;

  return 0;
}

