#!.venv/bin/python
"""
Execution for OSKAR.
"""

import sys
import os

sys.path.append("/scratch/mwaeor/ohrw/phd-programs/")
# pylint: disable=wrong-import-position
from oskareor.oskar_exec import BTAnalysisPipeline as btap
from oskareor.oskar_helpers import OSKARHelper as ohelp

SETONIX_DATA_DIR_SOFTWARE = "/software/projects/mwaeor/ohrw/"
SETONIX_DATA_DIR_SCRATCH = "/scratch/mwaeor/ohrw/"

# Execution Stage

MODEL = os.environ["OSKAR_MODEL"]
H5_LOCATION = os.environ["OSKAR_H5_LOCATION"]
DATA_DIR = SETONIX_DATA_DIR_SCRATCH + "/oskareor.data/oskar.data/" + "/" + MODEL + "/"

btap.run_oskar_on_model_timed(
    file=H5_LOCATION,
    oskar_parent_dir=SETONIX_DATA_DIR_SCRATCH,
    outpath=(
        DATA_DIR + MODEL + ohelp.TEMPLATE_FILE_TYPE_EXTENSIONS["ms"],
        DATA_DIR + MODEL + ohelp.TEMPLATE_FILE_TYPE_EXTENSIONS["vis"],
        DATA_DIR + MODEL + ohelp.TEMPLATE_FILE_TYPE_EXTENSIONS["fits"],
        DATA_DIR + MODEL + ohelp.TEMPLATE_FILE_TYPE_EXTENSIONS["uvfits"],
    ),
    oskar_mode="binary",
    oskar_exec=SETONIX_DATA_DIR_SCRATCH + ohelp.OSKAR_BIN,
    use_imager=True,
    load_osm=False,
    oskar_parent_dir=SETONIX_DATA_DIR_SCRATCH,
    convert_uvfits=True,
)
