#!/bin/sh

epsic -N 16 -s 1.1,0,0,1 -H 3
map2tga healpix.fits single+V.tga
rm healpix.fits

epsic -N 16 -s 1.1,0,0,-1 -H 3
map2tga healpix.fits single-V.tga
rm healpix.fits


epsic -N 16 -s 1.1,0,1,0 -H 3
map2tga healpix.fits single+U.tga
rm healpix.fits

epsic -N 16 -s 1.1,0,-1,0 -H 3
map2tga healpix.fits single-U.tga
rm healpix.fits


epsic -N 16 -s 1.1,1,0,0 -H 3
map2tga healpix.fits single+Q.tga
rm healpix.fits

epsic -N 16 -s 1.1,-1,0,0 -H 3
map2tga healpix.fits single-Q.tga
rm healpix.fits

