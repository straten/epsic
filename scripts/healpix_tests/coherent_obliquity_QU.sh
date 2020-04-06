#!/bin/bash

pi=3.14159265359

for ob in 000 020 040 060 080 100 120 140 160 180
do
  
  cos=`echo "c($ob*$pi/180)" | bc -l`
  sin=`echo "s($ob*$pi/180)" | bc -l`

  echo $ob $cos $sin
  epsic -N 16 -c 0.9 -s 2,$cos,$sin,0 -s B2,1,0,0 -H 3
  map2tga healpix.fits coherent_Q_U_${ob}.tga
  rm healpix.fits

done


