#!/bin/bash -l
#SBATCH --job-name=chips_ska
#SBATCH --output=chips.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=03:00:00
#SBATCH --partition=mwa-gpu
#SBATCH --account=mwaeor-gpu
#SBATCH --export=NONE
#SBATCH --gres=gpu:5
 
module use /software/projects/mwaeor/setonix/2024.05/modules/zen3/gcc/12.2.0
source /software/projects/mwaeor/ctrott/setonix/chips_2025_ska/env_variables_garrawarla.sh

export OUTPUTDIR='/scratch/mwaeor/ohrw/chips_run_stage'
export INPUTDIR='/scratch/mwaeor/ohrw/chips_run_stage'

cd $INPUTDIR

/software/projects/mwaeor/ctrott/setonix/chips_2025_ska/gridvisska /scratch/ohrw/uvfits_templates/uvfits_templates/slice_uvwplane.uvfits side

/software/projects/mwaeor/ctrott/setonix/chips_2025_ska/prepare_ska side 375 0 'yy' side 1

/software/projects/mwaeor/ctrott/setonix/chips_2025_ska/lssa_fg_ska side 375 50 'yy' 500. side 0 1

#fi