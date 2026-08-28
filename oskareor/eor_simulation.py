"""The EoR simulation model is for functions that interact with the 21cmFAST python wrapper."""

# System imports
import os

# Mathematics and calculations
import numpy as np
from astropy import units as u
from astropy.cosmology import z_at_value as getz

# Software
import py21cmfast as p21c  # type: ignore

# Local imports
from oskareor.skalow_calc import EoRCosmology as eorcosmo, SKAMath as omath, FileManager as ofmg, SKAString as ostr
from oskareor.oskar_helpers import OSKARHelper as ohelp


class Simulator:
    """A wrapper object for simulating the EoR with 21cmFast."""

    DEFAULT_SIMULATION_INPUTS = {
        "random_seed": 20250303,
        "HII_DIM": 512,
        "PERTURB_ON_HIGH_RES": True,
        "PHOTON_CONS_TYPE": "z-photoncons",
        "LOWRES_CELL_SIZE_MPC": 2,
        "N_THREADS": 8,
        "USE_CMB_HEATING": True,
        "USE_LYA_HEATING": True,
        "F_STAR10": -1.42,
        "ALPHA_STAR": 0.51,
        "F_ESC10": -1.08,
        "ALPHA_ESC": -0.46,
        "t_STAR": 0.34,
        "M_TURN": 8.46,
    }

    TEMPLATE_PRESETS = {
        "fiducial": (
            {"basic", "simple", "starter", "f", "s", "b"},
            "A simple fiducial model using a modified version of Park19 and the Gpc py21cmfast, templates.",
            {},
        )
    }

    @staticmethod
    def display_template_presets(print_presets: bool = True, filter_preset: str = "") -> dict:
        """Return and/or display all available templates, their names, and their descriptions. All templates except random, coeval1, and coeval2 are identical in all t-dimension voxels.

        :param print_presets: If true print the templates to console and return the dictionary. If false only return the dictionary.
        :param filter_preset: If true print/return the information for only one specific template entry, according to the given string.

        :return template_options: The dictionary containing all available templates or a dictionary of the specific desired template.
        """

        return ohelp.display_options(
            Simulator.TEMPLATE_PRESETS,
            print_options=print_presets,
            selection=filter_preset,
        )

    # FIXME: template config
    def __init__(
        self,
        template,
        toml_file="",
        z_ref=5.5,
        n_steps=1024,
        step_size_mpc=2,
        *override_p21c_templates,
        **override_simulation_inputs,
    ):
        # Define simulation input parameters and save to .toml file
        self.inputs = p21c.InputParameters.from_template(*override_p21c_templates, **override_simulation_inputs)

        self.toml_file = toml_file
        if self.toml_file != "" or self.toml_file is not None:
            p21c.write_template(self.inputs, self.toml_file)

            print(f"Parameters saved to {self.toml_file}")
            print("Generating cosmology and cube dimension information ...")

        # Creating an array of Mpc values similar to how 21cmFast does it

        # `cosmo_params` is dynamically provided by py21cmFAST, so use `getattr()` here to
        # avoid static-type warnings while still resolving the cosmology object at runtime.
        cosmo = getattr(self.inputs.cosmo_params, "cosmo")

        dz_ref = cosmo.comoving_distance(z_ref).to_value(u.Mpc)
        inc = step_size_mpc
        total = n_steps * inc
        end = dz_ref + total - inc

        print("21cmFast box depth parameters:")
        print("Start (cMpc)   =", dz_ref)
        print("Step (cMpc)    =", inc)
        print("Length (cMpc)  =", total)
        print("End (cMpc)     =", end)
        print("End redshift z =", getz(cosmo.comoving_distance, end * u.Mpc))

        self.lc_dist = ostr.generate(start=dz_ref, step=inc, num=n_steps) * u.Mpc

        print("Lightcone Distances: ", self.lc_dist, self.lc_dist.shape)
        print("Creating the cache directory ...")

        # Create the lightcone
        self.lcn = p21c.RectilinearLightconer(lc_distances=self.lc_dist, quantities=("brightness_temp", "density"))

    def default(self):
        """Default Park19+Gpc model. Differences between Park19 Model and this one:

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
        return Simulator(["Park19", "large"], **self.DEFAULT_SIMULATION_INPUTS)

    def run(self, temp_dir, output_file, backup_dir=""):
        """Run the 21cmFast lightcone simulation. Requires a lot of memory and time.

        :param temp_dir: The operating directory to save temporary simulation data.
        :param output_file: The h5 file to output simulation data.
        :param backup_dir: If not empty, the directory to save the temporary directory to.
        """
        cache = p21c.OutputCache(temp_dir)

        print("Running 21cmFast ...")

        # Run the lightcone
        lightcone = p21c.run_lightcone(lightconer=self.lcn, inputs=self.inputs, cache=cache, progressbar=True)

        print("Writing h5 data to file ...")

        # Save the data to both the temp folder and the phd_programs dir
        backup_file = cache.direc / "simulation_lightcone.h5"
        output_file = output_file + "simulation_lightcone.h5"
        lightcone.save(backup_file, clobber=True)
        lightcone.save(output_file, clobber=True)

        print("Saved lightcone to:", output_file, "; and ", backup_file)
