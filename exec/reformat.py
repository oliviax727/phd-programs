#!.venv/bin/python
"""
Run the regrid/reformatter/transformation program.
"""

import sys

sys.path.append("/home/olivia/Desktop/Programs/phd-programs/")
from oskareor.reformatter import SimulationReformatter as simref  # pylint: disable=wrong-import-position

SETONIX_DATA_DIR_SOFTWARE = "/software/projects/mwaeor/ohrw/"
SETONIX_DATA_DIR_SCRATCH = "/scratch/mwaeor/ohrw/"

# Reformatting

simref.generate_osm_from_h5(
    "/home/olivia/oskareor.data/simulations/project1/fiducial_lightcone_simulation.h5",
    osm_output="/home/olivia/oskareor.data/oskar.data/fiducial/fiducial_sky_model.osm",
    save_dynamic_settings="/home/olivia/oskareor.data/oskar.data/fiducial/fiducial_general_settings.ini",
)
