#! /bin/tcsh 

set SRCDIR=$HOME/Data/variance_of_varI

set delta = DELTA
set intensity = INTENSITY
set npoint = 100

set n_on  = `echo "$npoint * $delta" | bc -l`
set n_off = `echo "$npoint * (1.0-$delta)" | bc -l`

set I_on = `echo "1.0 + $intensity" | bc -l`
set I_off = "1.0"

cd CURDIR

set count=0

echo "on_mean on_var off_mean off_var src_mean src_var" > results.txt

while ( $count < 100 )

    epsic -s ${I_on},0,0,0 -N $n_on -d >& on.txt
    set on_var = `awk '$1=="var[0]" {print $3}' on.txt`
    set on_mean = `awk '$1=="mean[0]" {print $3}' on.txt`

    epsic -s ${I_off},0,0,0 -N $n_off -d >& off.txt
    set off_var = `awk '$1=="var[0]" {print $3}' off.txt`
    set off_mean = `awk '$1=="mean[0]" {print $3}' off.txt`

    set src_mean = `echo "$on_mean - $off_mean" | bc -l`
    set src_var = `echo "$on_var - $off_var - $src_mean * $off_mean" | bc -l`

    echo "$on_mean $on_var $off_mean $off_var $src_mean $src_var" >> results.txt

  @ count = $count + 1

end

