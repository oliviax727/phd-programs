#!.venv/bin/python
"""
Testing stage for OSKAR Regrid.
"""

# pylint: disable=line-too-long, unused-import, trailing-newlines
from matplotlib import pyplot as plt
from oskareor.oskar_exec import LoadDefaults, BTAnalysisPipeline
from oskareor.oskar_helpers import OSKARHelper as ohelp

# Testing stage

TEMPLATE = "flat"

BTAnalysisPipeline.run_oskar_on_model(
    template_preset=TEMPLATE,
    outpath=(
        ohelp.default_template_path(template_preset=TEMPLATE, file_type="ms"),
        ohelp.default_template_path(template_preset=TEMPLATE, file_type="vis"),
        ohelp.default_template_path(template_preset=TEMPLATE, file_type="fits"),
        ohelp.default_template_path(template_preset=TEMPLATE, file_type="uvfits"),
    ),
    oskar_mode="binary",
    oskar_exec="~" + ohelp.OSKAR_BIN,
)
