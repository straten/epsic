#!/bin/csh

g++ -o mixing `psrchive --cflags --ldflags --pgplot-libs` mixing.C
./mixing -u -D usb.eps/PS
./mixing -D dsb.eps/PS

