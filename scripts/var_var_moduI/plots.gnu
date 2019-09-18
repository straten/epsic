#!/usr/bin/env gnuplot

set terminal postscript enhanced monochrome solid "Times-Roman" 24

set output "stddev_of_sample_mean.eps"

# set logscale x
# set logscale y
# set xrange [3:6000]
unset key

set xlabel 'Predicted {/Symbol s}_S x 10^3 (Eqn A4)'
set ylabel 'Measured {/Symbol s}_S x 10^3' rotate by 90

f(x)=x
plot "sigma_mean_source.txt" using ($7*1000):($6*1000), f(x)
unset output

set output "stddev_of_sample_variance_of_sample_mean.eps"
set xlabel 'Predicted {/Symbol s}_v x 10^3'
set ylabel 'Measured {/Symbol s}_v x 10^3' rotate by 90
plot [0:0.42] "sigma_var_source.txt" using ($7*1000):($6*1000), f(x)
unset output

