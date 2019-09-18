#! /bin/tcsh 

set delta = DELTA
set intensity = INTENSITY
set M = SAMPLE
set beta = 0.9

set npoint = `echo "100 / $M" | bc -l`

set n_on  = `echo "$npoint * $delta" | bc -l`
set n_off = `echo "$npoint * (1.0-$delta)" | bc -l`

set I_on = `echo "1.0 + $intensity" | bc -l`
set I_off = "1.0"

set count=0

echo "on_mean on_var off_mean off_var src_mean src_var" > results.txt

while ( $count < 100 )

    epsic -S -s A:${intensity},0,0,0 -r A:$M -l A:$beta -s B:${I_off},0,0,0 -N $n_on -n $M -d >& on.txt
    set on_var = `awk '$1=="var[0]" {print $3}' on.txt`
    set on_mean = `awk '$1=="mean[0]" {print $3}' on.txt`

    epsic -s ${I_off},0,0,0 -N $n_off -n $M -d >& off.txt
    set off_var = `awk '$1=="var[0]" {print $3}' off.txt`
    set off_mean = `awk '$1=="mean[0]" {print $3}' off.txt`

    set src_mean = `echo "$on_mean - $off_mean" | bc -l`
    set src_var = `echo "$on_var - $off_var - ( $src_mean * $off_mean ) / $M" | bc -l`

    echo "$on_mean $on_var $off_mean $off_var $src_mean $src_var" >> results.txt

  @ count = $count + 1

end

