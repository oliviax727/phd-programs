#!.venv/bin/python
"""
Testing stage for OSKAR Regrid.
"""

# pylint: disable=line-too-long, unused-import
from matplotlib import pyplot as plt
from oskar_modules.oskar_exec import LoadDefaults, BTAnalysisPipeline

# Testing stage

#BTAnalysisPipeline.h5_box_to_datacube(None, template_preset="gaussian")

#BTAnalysisPipeline.h5_box_to_datacube("./regrid/osm_output/yuxiang1_zenith_osm", oskar_exec=RegridHelper.OSKAR_BIN, load_osm=True, oskar_mode="binary", oskar_telescope_model=RegridHelper.TELESCOPE)
