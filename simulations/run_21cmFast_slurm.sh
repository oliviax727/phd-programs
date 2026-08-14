#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=700GB
#SBATCH --time=12:00:00
#SBATCH --job-name=gpc-py21c-olivia
#SBATCH --account=oz113
#SBATCH --output=/fred/oz113/owalters/slurm.out

source /home/owalters/module-load.sh

export GSL_LIB="/fred/oz113/owalters/gsl-2.8/lib"
export GSL_INC="/fred/oz113/owalters/gsl-2.8/include"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/fred/oz113/owalters/gsl-2.8/lib"

cd /fred/oz113/owalters/phd-programs

./simulations/21cmfast_exec.py