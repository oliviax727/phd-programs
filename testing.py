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

flat_fields = [
    {"d": (400, 400, 1), "name": "flat400x2"},
    {"d": (500, 500, 1), "name": "flat500x2"},
    {"d": (600, 600, 1), "name": "flat600x2"},
]

for field in flat_fields:
    values = simref.mock_values(preset="flat400", d=field["d"])

    osm_output = "./oskar-chips/primary-beam-testing/" + field["name"] + ".osm"
    fits_output = "./oskar-chips/primary-beam-testing/" + field["name"] + ".fits"

    dynamic_settings = simref.generate_osm_from_simulation(values=values, z_ref=12, v=(2, 2, 2), osm_output=osm_output)

    btap.run_oskar_on_model(
        file=osm_output, outpath=("", "", fits_output, ""), load_osm=True, convert_uvfits=False, box_dim=field["d"]
    )
