#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=750GB
#SBATCH --time=4-00:00:00
#SBATCH --job-name=owalters-p21c
#SBATCH --account=oz113
#SBATCH --output=/fred/oz113/owalters/slurm.out
#SBATCH --tmp=1TB

echo "Initialising ..."

# Create Temporary Directory
run_date=$(date +"%F") || exit
export P21C_TEMPLATE="$1"

export P21C_TEMP_DIR="${JOBFS}/p21c-olivia-${P21C_TEMPLATE}_${run_date}"
export P21C_OUT_DIR="/fred/oz113/owalters"

mkdir -p "${P21C_TEMP_DIR}"

# TODO: Rsync command here!

echo "Writing cache data to: ${P21C_TEMP_DIR}"
echo "Loading modules ..."

# shellcheck disable=SC1091
source /home/owalters/module-load.sh

export GSL_LIB="/fred/oz113/owalters/gsl-2.8/lib"
export GSL_INC="/fred/oz113/owalters/gsl-2.8/include"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/fred/oz113/owalters/gsl-2.8/lib"

cd /fred/oz113/owalters/phd-programs || exit

echo "Running the simulation ..."

./exec/21cmfast.py

echo "Simulation complete, migrating data ..."

cp -r "$P21C_TEMP_DIR" "/fred/oz113/owalters/p21c.data/"

echo "Data transfer complete!"
