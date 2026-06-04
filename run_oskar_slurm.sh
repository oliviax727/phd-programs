#!/bin/bash -l
#SBATCH --job-name=oskareor-test
#SBATCH --output=/software/projects/mwaeor/ohrw
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=06:00:00
#SBATCH --partition=mwa-gpu
#SBATCH --account=mwaeor-gpu
#SBATCH --export=NONE
#SBATCH --gres=gpu:6

module load gcc-native/14.2
module load rust/1.85.0
module load hdf5/1.14.5-parallel-api-v112
module load casacore/3.5.0-wybngs5
module load cfitsio/4.4.0
module load rocm/6.4.1
module load libfabric/1.22.0
module load python/3.11.6

cd "/software/projects/mwaeor/ohrw/phd-programs"
.venv/bin/activate
./testing.py