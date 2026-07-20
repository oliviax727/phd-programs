#!/bin/bash -l
#SBATCH --job-name=ohrw_chips_ska
#SBATCH --output=chips.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=04:00:00
#SBATCH --partition=mwa-gpu
#SBATCH --account=mwaeor-gpu
#SBATCH --export=NONE
#SBATCH --gres=gpu:8

module use /software/projects/mwaeor/setonix/2024.05/modules/zen3/gcc/12.2.0
# shellcheck disable=SC1091
source /software/projects/mwaeor/ctrott/setonix/chips_2025_ska/env_variables_garrawarla.sh
module load cfitsio/4.4.0

export OUTPUTDIR='/scratch/mwaeor/ohrw/chips_run_stage/'
export INPUTDIR='/scratch/mwaeor/ohrw/chips_run_stage/'

cd $INPUTDIR || exit

rm -rf ../chips_run_stage/*
: > ../chips.out

mkdir -p ../2dps_png_templates
rm -rf ../2dps_png_templates/*

mkdir -p ../2dps_npz_templates
rm -rf ../2dps_npz_templates/*

chips_files=(/scratch/mwaeor/ohrw/uvfits_templates/uvfits_templates/*)

for template_file in "${chips_files[@]}"; do

    template_name="$(awk -v s="${template_file}" 'BEGIN { start = 39; end = index(s, "_uvw_plane.uvfits") - start; print substr(s, start, end) }')"

    /software/projects/mwaeor/ctrott/setonix/chips_2025_ska/gridvisska /scratch/mwaeor/ohrw/uvfits_templates/column_uvw_plane.uvfits side 15. 299792458.

    /software/projects/mwaeor/ctrott/setonix/chips_2025_ska/prepare_ska side 100 0 'yy' side 1

    /software/projects/mwaeor/ctrott/setonix/chips_2025_ska/lssa_fg_ska side 100 50 'yy' 500. chips.out 0 1

    # shellcheck disable=SC1091
    source /software/projects/mwaeor/ohrw/.venv/bin/activate

    plotchips_all.py --basedir . --plot_type 2D --chips_tag chips.out --polarisation YY --N_chan 100 --N_kperp 50 --max_power 1e10 --min_power 1e3 >> chips.log

    deactivate

    cp ./chips2D_yy_chips.out_crosspower.png ../2dps_png_templates/"${template_name}".png
    cp ./2D_coords_and_power.npz ../2dps_npz_templates/"${template_name}".npz

done