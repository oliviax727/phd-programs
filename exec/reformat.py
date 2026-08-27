#!.venv/bin/python
"""
Run the regrid/reformatter/transformation program.
"""

import sys

sys.path.append("/scratch/mwaeor/ohrw/phd-programs/")
from oskareor.oskar_exec import LoadDefaults as ldd  # pylint: disable=wrong-import-position

SETONIX_DATA_DIR_SOFTWARE = "/software/projects/mwaeor/ohrw/"
SETONIX_DATA_DIR_SCRATCH = "/scratch/mwaeor/ohrw/"

# Reformatting

ldd.reload_template_sky_models(oskar_parent_dir=SETONIX_DATA_DIR_SCRATCH)
