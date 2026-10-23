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

#cd /scratch/mwaeor/ohrw/phd-programs || exit

module use /software/projects/mwaeor/setonix/2024.05/modules/zen3/gcc/12.2.0
# shellcheck disable=SC1091
source /software/projects/mwaeor/ctrott/setonix/chips_2025_ska/env_variables_garrawarla.sh
module load cfitsio/4.4.0

# Run python file in mode 1 if CHUNK_DATA_FLAG flag is set
if [[ -n "${CHUNK_DATA_FLAG}" ]]; then
    rm -rf /scratch/mwaeor/ohrw/chips_run_stage

    mkdir -p /scratch/mwaeor/ohrw/chips_run_stage
    mkdir -p /scratch/mwaeor/ohrw/oskareor.data/oskar.data/fiducial/ps2ds

    ./exec/chips.py 1
fi

chunks=(/scratch/mwaeor/ohrw/chips_run_stage/PSA/chunks/*)

nchan=64
kperp=64
umax=1000
c=299792458

for chunk_file in "${chunks[@]}"; do
    # Run CHIPS
    echo "Running CHIPS ..."
    /software/projects/mwaeor/ctrott/setonix/chips_2025_ska/gridvisska "${chunk_file}" side 15.0 $c
    /software/projects/mwaeor/ctrott/setonix/chips_2025_ska/prepare_ska side $nchan 0 'yy' side 1
    /software/projects/mwaeor/ctrott/setonix/chips_2025_ska/lssa_fg_ska side $nchan $kperp 'yy' $umax chips.out 0 1

    # Extract chunk file data
    echo "CHIPS Finished, extracting chunk data ..."
    chunk="$(basename -- "${chunk_file}")"

    IFS='_' read -ra DATA <<< "${chunk}"
    f_chn="${DATA[1]}"
    f_ref="${DATA[2]}"
    index="${DATA[3]}"

    # Run python file in mode 2
    ./exec/chips.py 2 "${f_chn}" "${f_ref}" "${index}"

    echo "Chunk data extracted for ${index}"
done
