#!/bin/bash

current_dir=$(pwd)
cd /software/projects/mwaeor/ohrw/phd-programs/oskar_run_stage

# Run OSKAR
function run-oskar() {
    salloc --account=mwaeor-gpu --partition=gpu-dev --time=01:00:00 --gres=gpu:1 srun oskar-run.sh
}

# Run WSClean
function run-wsclean() {
    module load wsclean/3.4-idg
    bash wsclean.txt > test-wsc.out 2>&1
}

# Run all in chain
run-oskar && run-wsclean && ./compile-results.sh && cp -rf test-out /scratch/mwaeor/ohrw && echo "Test Complete" || echo "Test Faliure"

cd $current_dir