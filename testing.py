#!.venv/bin/python
"""
Testing stage for OSKAR Regrid.
"""

# pylint: disable=line-too-long, unused-import
from matplotlib import pyplot as plt
from oskareor.oskar_exec import LoadDefaults

# Testing stage

LoadDefaults.reload_template_sky_models(oskar_parent_dir = "/scratch/mwaeor/ohrw", update_which_templates = { "point" , "column", "gaussian"})
#LoadDefaults.reload_all(oskar_parent_dir = "/software/projects/mwaeor/ohrw", update_which_templates = { "point" })

