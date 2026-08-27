#!/bin/bash -l
#SBATCH --job-name=ohrw-oskareor-reformat
#SBATCH --output=/scratch/mwaeor/ohrw/slurm.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=mwa
#SBATCH --account=mwaeor
#SBATCH --export=NONE
#SBATCH --exclusive

cd /scratch/mwaeor/ohrw/phd-programs || exit
# shellcheck disable=SC1091
source .venv/bin/activate

# shellcheck disable=SC1091
source /software/projects/mwaeor/ohrw/install-scripts/oskar-install.sh

srun -N 1 -n 1 --exclusive ./exec/reformat.py
