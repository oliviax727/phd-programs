#!/bin/bash -l
#SBATCH --job-name=oskareor-test
#SBATCH --output=/software/projects/mwaeor/ohrw/slurm.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --partition=mwa-gpu
#SBATCH --account=mwaeor-gpu
#SBATCH --export=NONE
#SBATCH --gres=gpu:8
#SBATCH --exclusive

cd /software/projects/mwaeor/ohrw/phd-programs
source .venv/bin/activate

source ../install-scripts/oskar-install.sh

srun -N 1 -n 1 --exclusive ./testing.py
