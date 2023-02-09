#!/usr/bin/env bash

#PBS -P 2021_PIP                     
#PBS -l select=1:ncpus=18:mem=95GB   
#PBS -l walltime=400:00:00


cd $PBS_O_WORKDIR/

module load R/3.6.0

# let the R parallel library know it has 16 cores available
export MC_CORES=18

# ensure external libraries called by threads don't utilise additional parallelism
export OMP_NUM_THREADS=1


Rscript  R_SCRIPT_FILE


