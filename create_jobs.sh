#!/usr/bin/env bash

#PBS -P 2021_PIP 

cd $PBS_O_WORKDIR/

if [ ! -d "JOBS" ]; then
    mkdir "JOBS"
else
    echo "ERROR: the directory JOBS is already present!"
    exit
fi

for R_SCRIPT_FILE in *.R; do
    cp "job_template.sh" "JOBS/${R_SCRIPT_FILE}.sh" 
    sed -i "s/R_SCRIPT_FILE/${R_SCRIPT_FILE}/" "JOBS/${R_SCRIPT_FILE}.sh"
done

