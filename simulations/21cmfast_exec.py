#!.venv/bin/python

import py21cmfast as p21c
from tempfile import mkdtemp

print(f' Using 21cmFAST version {p21c.__version__}')

scin_inputs = p21c.InputParameters.from_template(
    ["Park19", "large"],
    random_seed = 20250303,
    HII_DIM = 512,
    LOWRES_CELL_SIZE_MPC = 2,
    N_THREADS = 8,
    USE_CMB_HEATING = True,
    USE_LYA_HEATING = True,
    F_STAR10 = -1.42,
    ALPHA_STAR = 0.51,
    F_ESC10 = -1.08,
    ALPHA_ESC = -0.46,
    t_STAR = 0.34,
    M_TURN = 8.46
)

scin_toml_file = "/fred/oz113/owalters/phd-programs/simulations/scin-park19.toml"

p21c.write_template(scin_inputs, scin_toml_file)

cache = p21c.OutputCache(mkdtemp())

lcn = p21c.RectilinearLightconer.between_redshifts(
    min_redshift=5.5,
    max_redshift=12.0,
    quantities=("brightness_temp", "density"),
    resolution=scin_inputs.simulation_options.cell_size
)

lightcone = p21c.run_lightcone(
    lightconer=lcn,
    inputs=scin_inputs,
    cache=cache,
    progressbar=True
)

filename = lightcone.save(cache.direc / "lightcone.h5", clobber=True)
filename2 = lightcone.save("/fred/oz113/owalters/phd-programs/simulations/simulation_lightcone.h5", clobber=True)

print("Saved lightcone to: " + filename)