#!/bin/bash -l
#SBATCH --job-name=ohrw-oskareor-sim
#SBATCH --output=/scratch/mwaeor/ohrw/slurm.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --partition=gpu-highmem
#SBATCH --account=mwaeor-gpu
#SBATCH --export=NONE
#SBATCH --gres=gpu:8
#SBATCH --exclusive

cd /scratch/mwaeor/ohrw/phd-programs
source .venv/bin/activate

source /software/projects/mwaeor/ohrw/install-scripts/oskar-install.sh

srun -N 1 -n 1 --exclusive "./exec/oskar.py"
