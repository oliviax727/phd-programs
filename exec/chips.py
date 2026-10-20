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

# Reformatting

cmgr.run_chips_on_model(
    SETONIX_DATA_DIR_SCRATCH + "oskareor.data/oskar.data/fiducial/fiducial_uvw_plane.uvfits",
    SETONIX_DATA_DIR_SCRATCH + "chips_run_stage",
    SETONIX_DATA_DIR_SCRATCH + "oskareor.data/oskar.data/fiducial/ps2ds",
)
