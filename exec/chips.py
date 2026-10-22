#!.venv/bin/python
"""
Run CHIPS.
"""

import sys

# sys.path.append("/home/olivia/Desktop/Programs/phd-programs/")
sys.path.append("/scratch/mwaeor/ohrw/phd-programs/")
from oskareor.power_spectra_analysis import ChipsManager as cmgr

SETONIX_DATA_DIR_SOFTWARE = "/software/projects/mwaeor/ohrw/"
SETONIX_DATA_DIR_SCRATCH = "/scratch/mwaeor/ohrw/"

# Locations

uvfits_loc = SETONIX_DATA_DIR_SCRATCH + "oskareor.data/oskar.data/fiducial/fiducial_uvw_plane.uvfits"
temp_parent_dir = SETONIX_DATA_DIR_SCRATCH + "chips_run_stage"
outdir = SETONIX_DATA_DIR_SCRATCH + "oskareor.data/oskar.data/fiducial/ps2ds"

# Arguments

args = sys.argv

match int(args[0]):
    case 0:
        print("Standard option selected")
        cmgr.run_chips_on_model(uvfits_loc, temp_parent_dir, outdir)

    case 1:
        print("Chunking option selected")
        cmgr.extract_uvfits_chunks(uvfits_loc, temp_parent_dir, outdir)

    case 2:
        print("Extraction option selected")

        [f_chn, f_ref, index] = args[1:4]

        cmgr.extract_and_save_2dps(
            outdir + f"fiducial_ps2d_data_chunk_{int(index)}.npz",
            temp_parent_dir=temp_parent_dir + "/PSA/",
            f_ref=float(f_ref),
            f_chn=float(f_chn),
        )
