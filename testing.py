#!.venv/bin/python
"""
Testing stage for OSKAR Regrid.
"""

# pylint: disable=line-too-long, unused-import, trailing-newlines
from matplotlib import pyplot as plt
from oskareor.oskar_exec import LoadDefaults as ldd, BTAnalysisPipeline as btap
from oskareor.reformatter import SimulationReformatter as simref
from oskareor.oskar_helpers import OSKARHelper as ohelp

# Testing stage

coevals = {"coeval1", "coeval2", "flat400"}

ldd.reload_all(update_which_templates=ldd.TEMPLATES - coevals)
