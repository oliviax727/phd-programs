#!.venv/bin/python
"""
Testing stage for OSKAR Regrid.
"""

# pylint: disable=line-too-long, unused-import
from matplotlib import pyplot as plt
from oskareor.oskar_exec import LoadDefaults

# Testing stage

LoadDefaults.reload_template_sky_models()
LoadDefaults.reload_template_oskar_sims(update_which_templates = [ "slice", "column" ])
