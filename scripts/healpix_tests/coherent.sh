#!/bin/sh

epsic -N 16 -c 0.9 -s 1,0,0,1 -s B1,0,0,-1 -H 3
map2tga healpix.fits coherent_V.tga
rm healpix.fits

epsic -N 16 -c 0.9 -s 3,0,0,1 -s B1,0,0,-1 -H 3
map2tga healpix.fits coherent+V.tga
rm healpix.fits

epsic -N 16 -c 0.9 -s 1,0,0,1 -s B3,0,0,-1 -H 3
map2tga healpix.fits coherent-V.tga
rm healpix.fits

epsic -N 16 -c 0.9 -s 1,0,1,0 -s B1,0,-1,0 -H 3
map2tga healpix.fits coherent_U.tga
rm healpix.fits

epsic -N 16 -c 0.9 -s 3,0,1,0 -s B1,0,-1,0 -H 3
map2tga healpix.fits coherent+U.tga
rm healpix.fits

epsic -N 16 -c 0.9 -s 1,0,1,0 -s B3,0,-1,0 -H 3
map2tga healpix.fits coherent-U.tga
rm healpix.fits

epsic -N 16 -c 0.9 -s 1,1,0,0 -s B1,-1,0,0 -H 3
map2tga healpix.fits coherent_Q.tga
rm healpix.fits

epsic -N 16 -c 0.9 -s 3,1,0,0 -s B1,-1,0,0 -H 3
map2tga healpix.fits coherent+Q.tga
rm healpix.fits

epsic -N 16 -c 0.9 -s 1,1,0,0 -s B3,-1,0,0 -H 3
map2tga healpix.fits coherent-Q.tga
rm healpix.fits

