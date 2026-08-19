#!.venv/bin/python
"""
Testing stage for OSKAR Regrid.
"""

# pylint: disable=line-too-long, unused-import, trailing-newlines
from matplotlib import pyplot as plt
from oskareor.oskar_exec import LoadDefaults as ldd, BTAnalysisPipeline as btap
from oskareor.reformatter import SimulationReformatter as simref
from oskareor.oskar_helpers import OSKARHelper as ohelp

SETONIX_DATA_DIR_SOFTWARE = "/software/projects/mwaeor/ohrw/"
SETONIX_DATA_DIR_SCRATCH = "/scratch/mwaeor/ohrw/"

# Testing stage

# Flat 512 is priority one
ldd.reload_all(update_which_templates=["flat512"])

# Then do everything else
ldd.reload_all(update_which_templates=ldd.TEMPLATES - {"flat512"})
