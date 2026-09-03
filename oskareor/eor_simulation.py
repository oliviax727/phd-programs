"""The EoR simulation model is for functions that interact with the 21cmFAST python wrapper."""

# Mathematics and calculations

from astropy import units as u
from astropy.cosmology import z_at_value as getz

# Software
import py21cmfast as p21c  # type: ignore

# Local imports
from oskareor.skalow_calc import FileManager as ofmg, SKAString as ostr
from oskareor.oskar_helpers import OSKARHelper as ohelp


class Simulator:
    """A wrapper object for simulating the EoR with 21cmFast."""

    class SimulationOptionValue(ohelp.OptionDictValue):
        """Simulation option dict value."""

        inputs: dict[str, any]
        z_ref: float
        n_steps: int
        step_size_mpc: float
        p21c_templates: list[str]

    SimulationOptions = ohelp.OptionDict[SimulationOptionValue]

    DEFAULT_SIMULATION_INPUTS = {
        "random_seed": 20250303,
        "HII_DIM": 512,
        "PERTURB_ON_HIGH_RES": True,
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

    DEFAULT_LIGHTCONE_QUANTITES = (
        "brightness_temp",
        "density",
        "ionisation_rate_G12",
        "kinetic_temperature",
    )

    TEMPLATE_PRESETS: SimulationOptions = {
        "test": {
            "aliases": {"t", "test"},
            "description": (
                "A simple fiducial model using a modified version of Park19 " + "and the Gpc py21cmfast, templates."
            ),
            "inputs": DEFAULT_SIMULATION_INPUTS | {"HII_DIM": 128, "PERTURB_ON_HIGH_RES": False},
            "z_ref": 5.5,
            "n_steps": 128,
            "step_size_mpc": 1,
            "p21c_templates": ["Park19", "small"],
        },
        "fiducial": {
            "aliases": {"basic", "simple", "starter", "f", "s", "b"},
            "description": (
                "A simple fiducial model using a modified version of Park19 " + "and the Gpc py21cmfast, templates."
            ),
            "inputs": DEFAULT_SIMULATION_INPUTS,
            "z_ref": 5.5,
            "n_steps": 1024,
            "step_size_mpc": 2,
            "p21c_templates": ["Park19", "large"],
        },
        "photoncons": {
            "aliases": {"pc", "phcons"},
            "description": (
                "A simple fiducial model using a modified version of Park19 " + "and the Gpc py21cmfast, templates."
            ),
            "inputs": DEFAULT_SIMULATION_INPUTS | {"PHOTON_CONS_TYPE": "z-photoncons"},
            "z_ref": 5.5,
            "n_steps": 1024,
            "step_size_mpc": 2,
            "p21c_templates": ["Park19", "large"],
        },
    }

    @staticmethod
    def display_template_presets(print_presets: bool = True, filter_preset: str = "") -> SimulationOptions:
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

    def __init__(
        self,
        oskareor_template="fiducial",
        name_override="",
        p21c_templates=None,
        z_ref=None,
        n_steps=None,
        step_size_mpc=None,
        **override_simulation_inputs,
    ):
        """Create a Simulator wrapper object.

        :param oskareor_template: The named template to base simulation params on, defaults to "fiducial".
        :param name_override: If not empty, use as the name for file saving purposes, defaults to "".
        :param z_ref: Refrence redshift, defaults to 5.5.
        :param n_steps: Number of voxels in the line-of-sight dimension, defaults to 1024.
        :param step_size_mpc: Size of each line-of-sight dimension voxel, defaults to 2.
        """

        # Get oskareor template, default is fiducial
        self.template_name = oskareor_template
        self.template = self.display_template_presets(False, self.template_name)

        self.name = self.template_name if name_override == "" else name_override

        if z_ref is None:
            z_ref = self.template[self.template_name]["z_ref"]

        if n_steps is None:
            n_steps = self.template[self.template_name]["n_steps"]

        if step_size_mpc is None:
            step_size_mpc = self.template[self.template_name]["step_size_mpc"]

        if p21c_templates is None:
            p21c_templates = self.template[self.template_name]["p21c_templates"]

        print(p21c_templates, (self.template[self.template_name]["inputs"] | override_simulation_inputs))
        # Define simulation input parameters
        self.inputs = p21c.InputParameters.from_template(
            p21c_templates, **(self.template[self.template_name]["inputs"] | override_simulation_inputs)
        )

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

        # Create empty output data
        self.toml_file = self.file_name = self.output_file = self.backup_file = ""

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

    def run(self, temp_dir, out_dir, quantities=None):
        """Run the 21cmFast lightcone simulation. WARNING! Requires a lot of memory and time!

        :param temp_dir: The operating directory to save temporary simulation data.
        :param output_dir: The directory to output simulation data.
        :param quantities: The global quantities that should be simulated in-detail, by default this will be the `DEFAULT_LIGHTCONE_QUANTITES` object.
        :param backup_dir: If not empty, the directory to save the temporary directory to.
        """
        if quantities is None:
            quantities = self.DEFAULT_LIGHTCONE_QUANTITES

        cache = p21c.OutputCache(temp_dir)

        # Saving inputs to cache
        self.toml_file = ofmg.expand_path(temp_dir + "/" + self.name + "_input_params.toml")

        p21c.write_template(self.inputs, self.toml_file)

        print(f"Parameters saved to {self.toml_file}")

        print("Running 21cmFast ...")

        # Run the lightcone
        lightcone = p21c.run_lightcone(lightconer=self.lcn, inputs=self.inputs, cache=cache, progressbar=True)

        print("Writing h5 data to file ...")

        # Save the data to both the temp folder and the phd_programs dir
        self.file_name = self.name + "_lightcone_simulation.h5"
        self.backup_file = ofmg.expand_path(str(cache.direc) + "/" + self.file_name)
        lightcone.save(self.backup_file, clobber=True)
        self.output_file = ofmg.expand_path(out_dir + "/" + self.file_name)
        lightcone.save(self.output_file, clobber=True)

        print("Saved lightcone to:", self.output_file, "; and ", self.backup_file)
