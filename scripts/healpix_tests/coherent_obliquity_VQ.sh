#!/bin/bash

pi=3.14159265359

for ob in 000 020 040 060 080 100 120 140 160 180
do
  
  cos=`echo "c($ob*$pi/180)" | bc -l`
  sin=`echo "s($ob*$pi/180)" | bc -l`

  echo $ob $cos $sin
  epsic -N 16 -c 0.9 -s 2,$sin,0,$cos -s B2,0,0,-1 -H 3
  map2tga healpix.fits coherent_V_Q_${ob}.tga
  rm healpix.fits

done

exit -1


  echo $Q
  epsic -N 16 -c 0.9 -s 2,0,0,1 -s B2,$Q,0,-1 -H 3
  map2tga healpix.fits coherent_V_Q_${Q}.tga
  rm healpix.fits


