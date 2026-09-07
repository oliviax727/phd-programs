#!.venv/bin/python
"""Run 21cmFast."""

# Import simulation module
import sys
import os

sys.path.append("/fred/oz113/owalters/phd-programs/")
from oskareor.eor_simulation import Simulator  # pylint: disable=wrong-import-position

print("Getting environment variables ...")

TEMP_DIR = os.environ["P21C_TEMP_DIR"]
OUT_DIR = os.environ["P21C_OUT_DIR"]
TEMPLATE = os.environ["P21C_TEMPLATE"]

sim = Simulator(oskareor_template=TEMPLATE)

sim.run(TEMP_DIR, OUT_DIR)
