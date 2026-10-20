#!/bin/bash -l
#SBATCH --job-name=ohrw_chips_ska
#SBATCH --output=/scratch/mwaeor/ohrw/chips.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=12:00:00
#SBATCH --partition=mwa-gpu
#SBATCH --account=mwaeor-gpu
#SBATCH --export=NONE
#SBATCH --gres=gpu:8

cd /scratch/mwaeor/ohrw/phd-programs || exit

module use /software/projects/mwaeor/setonix/2024.05/modules/zen3/gcc/12.2.0
# shellcheck disable=SC1091
source /software/projects/mwaeor/ctrott/setonix/chips_2025_ska/env_variables_garrawarla.sh
module load cfitsio/4.4.0

mkdir -p /scratch/mwaeor/ohrw/chips_run_stage
mkdir -p /scratch/mwaeor/ohrw/oskareor.data/oskar.data/fiducial/ps2ds

./exec/chips.py
