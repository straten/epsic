#!/bin/bash

pi=3.14159265359

for ob in 000 020 040 060 080 100 120 140 160 180
do
  
  cos=`echo "c($ob*$pi/180)" | bc -l`
  sin=`echo "s($ob*$pi/180)" | bc -l`

  echo $ob $cos $sin
  epsic -N .001953125 -c 0.99 -s 2,$sin,0,$cos -s B2,0,0,-1 -f -n 1024
  mv stokes.txt ${ob}.txt

done


