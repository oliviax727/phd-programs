#!.venv/bin/python
"""Test run 21cmFast."""

# Import simulation module
import sys

# import os

sys.path.append("/home/olivia/Desktop/Programs/phd-programs/")
from oskareor.eor_simulation import Simulator  # pylint: disable=wrong-import-position

print("Getting environment variables ...")

for template in ["q25-nospin"]:
    print("Simulating template:", template)

    sim = Simulator(oskareor_template=template)

    sim.run("/tmp/", "/home/olivia/oskareor.data/simulations/project1/")

    del sim
