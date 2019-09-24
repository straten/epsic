#! /bin/csh

# the duty cycle, delta
foreach delta (0.05 0.10 0.20 0.40)

 # the intensity ratio, r
 foreach intensity (0.001 0.01 0.1 1.0)

  # the sub-sample size, m
  foreach sample ( 10 20 50 )

   echo "********************* " delta: $delta intensity: $intensity sample: $sample

   # index starts at ten simply so that all directories have double-digit numbers in names
   set index = 10

   # divide job into 32 parts
   while ( $index < 42 )

    set DIR=results/data_${index}_${delta}_${intensity}_${sample}
    mkdir -p $DIR

    sed -e "s|DELTA|$delta|g" -e "s|INTENSITY|$intensity|g" -e "s|SAMPLE|$sample|g" single.csh > $DIR/trial.csh

    cd $DIR
    sbatch --nodes=1 --time=60 trial.csh
    cd -

    @ index = $index + 1

   end

  end 

 end

end

