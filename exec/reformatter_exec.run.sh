#!/bin/bash -l
#SBATCH --job-name=oskareor-test
#SBATCH --output=/scratch/mwaeor/ohrw/slurm.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=mwa-gpu
#SBATCH --account=mwaeor-gpu
#SBATCH --export=NONE
#SBATCH --gres=gpu:12
#SBATCH --exclusive

cd /scratch/mwaeor/ohrw/phd-programs
source .venv/bin/activate

source /software/projects/mwaeor/ohrw/install-scripts/oskar-install.sh

srun -N 1 -n 1 --exclusive ./testing.py
