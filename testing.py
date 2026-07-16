#!.venv/bin/python
"""
Testing stage for OSKAR Regrid.
"""

# pylint: disable=line-too-long, unused-import, trailing-newlines
from matplotlib import pyplot as plt
from oskareor.oskar_exec import LoadDefaults

# Testing stage

coevals = { "coeval1", "coeval2" }
templates = LoadDefaults.TEMPLATES - coevals
LoadDefaults.reload_all(update_which_templates=templates)
LoadDefaults.reload_all(update_which_templates=coevals)
