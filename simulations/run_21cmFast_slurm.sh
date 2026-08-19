#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=300GB
#SBATCH --time=24:00:00
#SBATCH --job-name=gpc-py21c-olivia
#SBATCH --account=oz113
#SBATCH --output=/fred/oz113/owalters/slurm.out
#SBATCH --tmp=2TB

# Create Temporary Directory
export P21C_TEMP_DIR=$JOBFS/p21c

# shellcheck disable=SC1091
source /home/owalters/module-load.sh

export GSL_LIB="/fred/oz113/owalters/gsl-2.8/lib"
export GSL_INC="/fred/oz113/owalters/gsl-2.8/include"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/fred/oz113/owalters/gsl-2.8/lib"

cd /fred/oz113/owalters/phd-programs || exit

./simulations/21cmfast_exec.py

tar -czf /fred/oz113/owalters/p21c.data.tar.gz "$P21C_TEMP_DIR"
