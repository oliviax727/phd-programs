#!.venv/bin/python
"""
Run the regrid/reformatter/transformation program.
"""

import os
import sys

sys.path.append(os.getcwd())
from oskareor.reformatter import SimulationReformatter as simref  # pylint: disable=wrong-import-position

SETONIX_DATA_DIR_SOFTWARE = "/software/projects/mwaeor/ohrw/"
SETONIX_DATA_DIR_SCRATCH = "/scratch/mwaeor/ohrw/"

parent = SETONIX_DATA_DIR_SCRATCH  # "~/"  # SETONIX_DATA_DIR_SCRATCH

# Reformatting

simref.generate_osm_from_h5(
    parent + "oskareor.data/simulations/project1/fiducial_lightcone_simulation.h5",
    osm_output=parent + "oskareor.data/oskar.data/fiducial/fiducial_sky_model.osm",
    save_dynamic_settings=parent + "oskareor.data/oskar.data/fiducial/fiducial_general_settings.ini",
)
