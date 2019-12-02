#! /bin/csh

cd results

rm -f source_results.txt off_results.txt product_results.txt

foreach delta (0.05 0.10 0.20 0.40)

 foreach intensity (0.001 0.01 0.1 1.0)

  foreach sample ( 10 20 50 )

    set file=results_${delta}_${intensity}_${sample}.txt
    echo creating results/$file

    cat data_*_${delta}_${intensity}_${sample}/results.txt | grep -v on > $file

  end
 end
end

