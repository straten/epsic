#! /bin/csh

set SRCDIR=$HOME/Data/variance_of_var_meanI
set OUTDIR=/fred/oz002/wvanstra/variance_of_var_meanI

rm -f slurm-*.out

foreach delta (0.05 0.10 0.20 0.40)
	
  foreach intensity (0.001 0.01 0.1 1.0)

    foreach sample ( 10 20 50 )

    echo "********************* " delta: $delta intensity: $intensity sample: $sample

  set index = 10

  while ( $index < 42 )

    set CURDIR=$OUTDIR/data_${index}_${delta}_${intensity}_${sample}
    mkdir -p $CURDIR

    sed -e "s|DELTA|$delta|g" -e "s|INTENSITY|$intensity|g" -e "s|SAMPLE|$sample|g" -e "s|CURDIR|$CURDIR|g" single.csh > $CURDIR/trial.csh

    sbatch --nodes=1 --time=60 $CURDIR/trial.csh

    @ index = $index + 1

  end

    end 

  end

end

