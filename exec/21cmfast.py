#!.venv/bin/python
"""Differences between Park19 Model and this one:

HII_DIM=512
LOWRES_CELL_SIZE_MPC=2
PERTURB_ON_HIGH_RES=True
PHOTON_CONS_TYPE=z-photoncons

From Meiksin 2021 and Reis et al. 2021:
- USE_CMB_HEATING = True
- USE_LYA_HEATING = True

From Qin et al. 2024 (Table 1):

- F_STAR10 = -1.42
- ALPHA_STAR = 0.51
- F_ESC10 = -1.08
- ALPHA_ESC = -0.46
- t_STAR = 0.34
- M_TURN = 8.46

This assumes a constant f_esc. Pursuing a non-constant f_esc will be optional later on.
"""

# Module imports and environment variables
import os

import numpy as np
import py21cmfast as p21c
from astropy import units as u
from astropy.cosmology import z_at_value as getz

print("Getting environment variables ...")

TEMP_DIR = os.environ["P21C_TEMP_DIR"]
OUT_DIR = os.environ["P21C_OUT_DIR"]

print(f" Using 21cmFAST version {p21c.__version__}")
print("Generating inputs object ...")

# Define simulation input parameters and save to .toml file
inputs = p21c.InputParameters.from_template(
    ["Park19", "large"],
    random_seed=20250303,
    HII_DIM=512,
    PERTURB_ON_HIGH_RES=True,
    PHOTON_CONS_TYPE="z-photoncons",
    LOWRES_CELL_SIZE_MPC=2,
    N_THREADS=8,
    USE_CMB_HEATING=True,
    USE_LYA_HEATING=True,
    F_STAR10=-1.42,
    ALPHA_STAR=0.51,
    F_ESC10=-1.08,
    ALPHA_ESC=-0.46,
    t_STAR=0.34,
    M_TURN=8.46,
)

PARAMS_TOML = "/fred/oz113/owalters/phd-programs/simulations/scin-park19.toml"

p21c.write_template(inputs, PARAMS_TOML)

print(f"Parameters saved to {PARAMS_TOML}")
print("Generating cosmology and cube dimension information ...")

# Creating an array of Mpc values similar to how 21cmFast does it

# `cosmo_params` is dynamically provided by py21cmFAST, so use `getattr()` here to
# avoid static-type warnings while still resolving the cosmology object at runtime.
cosmo = getattr(inputs.cosmo_params, "cosmo")
Z_REF = 5.5
N_STEPS = 1024
STEP_MPC = 2 * u.Mpc

dz_ref = cosmo.comoving_distance(Z_REF).to_value(u.Mpc)
inc = STEP_MPC.to_value(u.Mpc)
total = N_STEPS * inc
end = dz_ref + total - inc

print("21cmFast box depth parameters:")
print("Start (cMpc)   =", dz_ref)
print("Step (cMpc)    =", inc)
print("Length (cMpc)  =", total)
print("End (cMpc)     =", end)
print("End redshift z =", getz(cosmo.comoving_distance, end * u.Mpc))


def generate(start=0, step=1, num=1):
    """Create an array of linearly-increasing consecutive entries with a defined start, step size, and step count. Similar to np.linspace and np.arrange."""
    return start + np.arange(0, num) * step


lc_dist = generate(start=dz_ref, step=inc, num=N_STEPS) * u.Mpc

print("Lightcone Distances: ", lc_dist, lc_dist.shape)
print("Creating the cache directory ...")

# Create the lightcone and cache
cache = p21c.OutputCache(TEMP_DIR)
lcn = p21c.RectilinearLightconer(lc_distances=lc_dist, quantities=("brightness_temp", "density"))

print("Running 21cmFast ...")

# Run the lightcone
lightcone = p21c.run_lightcone(lightconer=lcn, inputs=inputs, cache=cache, progressbar=True)

print("Writing h5 data to file ...")

# Save the data to both the temp folder and the phd_programs dir
backup_file = cache.direc / "simulation_lightcone.h5"
output_file = OUT_DIR + "simulation_lightcone.h5"
lightcone.save(backup_file, clobber=True)
lightcone.save(output_file, clobber=True)

print("Saved lightcone to:", output_file, "; and ", backup_file)
