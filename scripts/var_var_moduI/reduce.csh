#! /bin/csh

cd results

rm -f sigma_*.txt covar_*.txt

set N = 104857600
set beta = 0.9

echo "delta r m N <Up> sigUp predict approx" | tee sigma_var_source.txt

foreach delta (0.05 0.10 0.20 0.40)

  foreach intensity (0.001 0.01 0.1 1.0)

   foreach sample ( 10 20 50 )

    set file=results_${delta}_${intensity}_${sample}.txt

    ## Compute the standard deviation of the second noncentral moment of the off-pulse phase bins
    set mean = `awk '{count++; mu2=$4+$3*$3; sum+=mu2; sumsq+=mu2*mu2} END{mu=sum/count; var=sumsq/count-mu*mu; print count, mu, sqrt(var)}' $file`

    set sigma_mu2_off = `echo "sqrt( (2.0 + 1.0/(2.0 * $sample)) / ( $N * (1.0-$delta) ) )" | bc -l`

    echo ${delta} ${intensity} ${sample} ${mean} $sigma_mu2_off | tee >> sigma_mu2_off.txt

    
    ## Compute the standard deviation of the second noncentral moment of the on-pulse phase bins
    set mean = `awk '{count++; mu2=$2+$1*$1; sum+=mu2; sumsq+=mu2*mu2} END{mu=sum/count; var=sumsq/count-mu*mu; print count, mu, sqrt(var)}' $file`

    set var_on = `echo "$beta * $intensity" | bc -l`
    set sigma_mu2_on = `echo "$var_on^2 + sqrt( (2.0 + 1.0/(2.0 * $sample)) / ( $N * $delta ) )" | bc -l`

    echo ${delta} ${intensity} ${sample} ${mean} $sigma_mu2_on | tee >> sigma_mu2_on.txt


    # (Square root of) co-variance between first and second noncentral moments of off-pulse intensities
    set mean = `awk '{count++; p=$4+$3*$3; sump+=p; sumpsq+=p*p; sumI+=$3; sumIsq+=$3*$3; cov+=p*$3} END{mup=sump/count; varp=sumpsq/count-mup*mup; muI=sumI/count; varI=sumIsq/count-muI*muI; rho=cov/count-mup*muI; print sqrt(rho)}' $file`

    set covar_off = `echo "sqrt( 1.0 / ($N * (1.0-$delta)) )" | bc -l`
    echo ${delta} ${intensity} ${sample} ${mean} $covar_off | tee >> covar_off.txt


    # (Square root of) co-variance between first and second noncentral moments of on-pulse intensities
    set mean = `awk '{count++; p=$2+$1*$1; sump+=p; sumpsq+=p*p; sumI+=$1; sumIsq+=$1*$1; cov+=p*$1} END{mup=sump/count; varp=sumpsq/count-mup*mup; muI=sumI/count; varI=sumIsq/count-muI*muI; rho=cov/count-mup*muI; print sqrt(rho)}' $file`

    set I_on = `echo 1.0+$intensity | bc -l`
    set covar_on = `echo "sqrt( $I_on^3 / ($N * $delta) )" | bc -l`
    echo ${delta} ${intensity} ${sample} ${mean} $covar_on | tee >> covar_on.txt


    # Compute the standard deviation of the source mean intensity estimator
    set mean = `awk '{count++; p=$1-$3; sum+=p; sumsq+=p*p} END{mu=sum/count; var=sumsq/count-mu*mu; print count, mu, sqrt(var)}' $file`

    set I_S = $intensity
    set r = $I_S
    set sigma_mean_source = `echo "sqrt( (1 + 2*(1.0-$delta) * $r * (1 + $r/2) ) / ( 2.0 * $N * $delta * (1.0-$delta)) )" | bc -l`

    echo ${delta} ${intensity} ${sample} ${mean} $sigma_mean_source | tee >> sigma_mean_source.txt


    ## Compute the standard deviation of the source variance estimator

    set d_don  = `echo "-2.0 * $I_on - 1.0/$sample" | bc -l`
    set d_doff = `echo "2.0 - ($I_on - 2.0)/$sample" | bc -l`

    set var_on = `echo "${I_on}^2 / ( 2.0 * $N * $delta )" | bc -l`
    set var_off = `echo "1.0 / ( 2.0 * $N * (1.0-$delta) )" | bc -l`

    set term1 = `echo "$d_don^2 * $var_on" | bc -l`
    set term2 = `echo "$d_doff^2 * $var_off" | bc -l`
    set term3 = `echo "2.0 * $d_don * $covar_on^2" | bc -l`
    set term4 = `echo "-2.0 * $d_doff * $covar_off^2" | bc -l`

    ## Compute the standard deviation of the source variance estimator

    set mean = `awk -vm=$sample '{count++; up=$2-$4-($5*$3)/m; sum+=up; sumsq+=up*up} END{mu=sum/count; var=sumsq/count-mu*mu; print count, mu, sqrt(var)}' $file`

    set sigma_var_source = `echo "sqrt( ${sigma_mu2_on}^2 + ${sigma_mu2_off}^2 + $term1 + $term2 + $term3 + $term4 )" | bc -l`

    set approximation = `echo "sqrt( (1.0 + 1.0/$sample) / ( 2.0 * $N * $delta * (1.0-$delta) ) )" | bc -l`

    echo ${delta} ${intensity} ${sample} ${mean} $sigma_var_source $approximation | tee -a sigma_var_source.txt

  end
 end
end

