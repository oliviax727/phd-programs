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

make_chunks=0
is_template=0
preset_name="fiducial"

# Extract arguments
while [ $# -gt 0 ]; do
    case $1 in
        -c)
        make_chunks=0
        ;;
        -t)
        is_template=1
        ;;
        *)
        preset_name="$1"
        ;;
    esac
    shift
done

# Make the output directory
if [[ "${is_template}" == "1" ]]; then
    mkdir -p "/scratch/mwaeor/ohrw/oskareor.data/2dps_npz_templates/${preset_name}_ps2ds/"
elif [[ "${is_template}" == "0" ]]; then
    mkdir -p "/scratch/mwaeor/ohrw/oskareor.data/oskar.data/${preset_name}/ps2ds/"
fi

# Run python file in mode 1 or 3 if the chunk option is set
if [[ "${make_chunks}" == "0" ]]; then
    ./exec/chips.py "$((is_template * 2 + 1))" "${preset_name}"
fi

chunks=(/scratch/mwaeor/ohrw/PSA/chunks/*)

nchan=64
kperp=64
umax=1000
c=299792458

export OUTPUTDIR='/scratch/mwaeor/ohrw/PSA/chips_output/'
export INPUTDIR='/scratch/mwaeor/ohrw/PSA/chips_output/'

for chunk_file in "${chunks[@]}"; do
    # Run CHIPS
    echo "Running CHIPS on $chunk_file"
    cd $INPUTDIR || exit

    /software/projects/mwaeor/ctrott/setonix/chips_2025_ska/gridvisska "${chunk_file}" side 15.0 $c
    /software/projects/mwaeor/ctrott/setonix/chips_2025_ska/prepare_ska side $nchan 0 'yy' side 1
    /software/projects/mwaeor/ctrott/setonix/chips_2025_ska/lssa_fg_ska side $nchan $kperp 'yy' $umax chips.out 0 1

    # Extract chunk file data
    echo "CHIPS Finished processing ${chunk_file}, extracting chunk data ..."
    chunk="$(basename -- "${chunk_file}")"

    IFS='_' read -ra DATA <<< "${chunk}"
    f_chn="${DATA[1]}"
    f_ref="${DATA[2]}"
    index="${DATA[3]}"

    cd /scratch/mwaeor/ohrw/phd-programs || exit

    # Run python file in mode 2 or 4
    ./exec/chips.py "$((is_template * 2 + 2))" "${preset_name}" "${f_chn}" "${f_ref}" "${index}"

    echo "Chunk data extracted for ${index}"
done

