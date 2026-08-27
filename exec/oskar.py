#!.venv/bin/python
"""
Execution for OSKAR.
"""

import sys
sys.path.append("/scratch/mwaeor/ohrw/phd-programs/")
from oskareor.oskar_exec import LoadDefaults as ldd

SETONIX_DATA_DIR_SOFTWARE = "/software/projects/mwaeor/ohrw/"
SETONIX_DATA_DIR_SCRATCH = "/scratch/mwaeor/ohrw/"

# Execution Stage

ldd.reload_template_oskar_sims(oskar_parent_dir=SETONIX_DATA_DIR_SCRATCH)
