#!/usr/bin/env bash

#PBS -P 2021_PIP
 
cd $PBS_O_WORKDIR/

for JOB_SCRIPT in JOBS/*
 do
     qsub -N ${JOB_SCRIPT}  ${JOB_SCRIPT}
done 

