#! /bin/csh

set OUTDIR=/fred/oz002/wvanstra/variance_of_varI

rm -f source_results.txt off_results.txt product_results.txt

set N = 104857600

foreach delta (0.05 0.10 0.20 0.40)

  foreach intensity (0.001 0.01 0.1 1.0)

    set file=results_${delta}_${intensity}.txt

    grep -v on_mean $OUTDIR/data_*_${delta}_${intensity}/results.txt > $file

    ## Variance of source variance estimator

    set mean = `awk '{count++; sum+=$6; sumsq+=$6*$6} END{mu=sum/count; var=sumsq/count-mu*mu; print count, mu, sqrt(var)}' $file`

    set prediction = `echo "sqrt( 7.0 / ( 4.0 * $N * $delta * (1.0-$delta) ) )" | bc -l`
    
    echo ${delta} ${intensity} ${mean} $prediction | tee >> source_results.txt

    ## Variance of off-pulse phase bins
set mean = `awk '{count++; sum+=$4; sumsq+=$4*$4} END{mu=sum/count; var=sumsq/count-mu*mu; print count, mu, sqrt(var)}' $file`

    set prediction = `echo "sqrt( 5.0 / ( 4.0 * $N * (1.0-$delta) ) )" | bc -l`

    echo ${delta} ${intensity} ${mean} $prediction | tee >> off_results.txt

    ## Variance of on-pulse phase bins
set mean = `awk '{count++; sum+=$2; sumsq+=$2*$2} END{mu=sum/count; var=sumsq/count-mu*mu; print count, mu, sqrt(var)}' $file`

    set I_on = `echo 1.0+$intensity | bc -l`
    set prediction = `echo "$I_on * $I_on * sqrt( 5.0 / ( 4.0 * $N * $delta ) )" | bc -l`

    echo ${delta} ${intensity} ${mean} $prediction | tee >> on_results.txt

    # product of intensities
    set mean = `awk '{count++; p=($1-$3)*$3; sum+=p; sumsq+=p*p} END{mu=sum/count; var=sumsq/count-mu*mu; print count, mu, sqrt(var)}' $file`

    set I_S = $intensity
    set r = $I_S
    set prediction = `echo "sqrt( $I_S * $I_S / ( 2.0 * $N * (1.0-$delta) ) + ( 1 + 2 * (1.0-$delta) * $r * ( 2 + 0.5*$r ) ) / ( 2.0 * $N * $delta ) )" | bc -l`

    echo ${delta} ${intensity} ${mean} $prediction | tee >> product_results.txt

  end
end

