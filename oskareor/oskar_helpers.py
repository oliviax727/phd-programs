"""The oskar_helpers module contains a series of helper functions and variables for simulation-to-power-spectra pipeline. Requires an existing and robust oskareor.data file in the home directory."""

# Mathematics and calculations
# Performance Handling
from timeit import default_timer as timer
from typing import Any, Callable

import h5py
import numpy as np

# In-Project Objects
from oskareor.skalow_calc import FileManager as ofmg
from oskareor.skalow_calc import SKAString as ostr


# OSKAR EOR helper class - dependant on ~/oskareor.data
class OSKARHelper:
    """Helper functions and constants specific to handling OSKAR and the Reformatter."""

    # Default paths - Primary
    OSKAR_SIF: str = "/oskareor.data/OSKAR-2.12.2-Python3.sif"
    OSKAR_BIN: str = "/oskareor.data/bin/"
    TELESCOPE: str = "/oskareor.data/SKA-Low_telescope_models/SKA-Low_AAstar_original_rigid-rotation.tm"

    # Default paths - Templates
    TEMPLATE_FILE_TYPE_EXTENSIONS: dict[str, str] = {
        "osm": "_sky_model.osm",
        "ini": "_general_settings.ini",
        "vis": "_visibilities.vis",
        "ms": "_measurement_set.ms",
        "fits": "_datacube.fits",
        "uvfits": "_uvw_plane.uvfits",
    }

    @staticmethod
    def default_template_path(template_preset: str, oskar_parent_dir: str = "~", file_type: str = "") -> str:
        """Constructs a default template path.

        :param preset: Mock brightness temperature array format. Run SimulationReformatter.display_template_presets for more information.
        :param oskar_parent_dir: The directory containing the oskareor.data folder (default is the home folder).
        :param file_type: The file type extension to construct the path around.

        :returns template_path: The path to a specified template file.
        """

        return (
            oskar_parent_dir
            + "/oskareor.data/"
            + file_type
            + "_templates/"
            + template_preset
            + OSKARHelper.TEMPLATE_FILE_TYPE_EXTENSIONS[file_type]
        )

    # Define default settings
    DEFAULT_INTERFEROMETER_SETTINGS: dict = {
        "General": {"app": "oskar_sim_interferometer"},
        "simulator": {"use_gpus": True, "write_status_to_log_file": True},
        "observation": {
            "num_time_steps": 24,
            # Below values are modified dynamically =>
            "start_frequency_hz": 200e6,
            "num_channels": 100,
            "frequency_inc_hz": 140e3,
            "phase_centre_ra_deg": 0.0,
            "phase_centre_dec_deg": -27.0,
            # <=
            "length": "4:00:00.00",
            "start_time_utc": "2025-03-03 03:30:00.00",
        },
        "telescope": {
            "input_directory": "BTA/telescope_model.tm",
            "aperture_array/element_pattern/enable_numerical": False,
            # "aperture_array/array_pattern/element/x_gain": 1.0,
            # "aperture_array/array_pattern/element/y_gain": 1.0,
            # "aperture_array/array_pattern/element/x_gain_error_time": 0.0015057,
            # "aperture_array/array_pattern/element/y_gain_error_time": 0.0015057,
            # "aperture_array/array_pattern/element/x_phase_error_fixed_deg": 0.0,
            # "aperture_array/array_pattern/element/y_phase_error_fixed_deg": 0.0,
            # "aperture_array/array_pattern/element/x_phase_error_time_deg": 0.0015057,
            # "aperture_array/array_pattern/element/y_phase_error_time_deg": 0.0015057
        },
        "interferometer": {
            "oskar_vis_filename": "BTA/oskar_output/vis.vis",
            "ms_filename": "BTA/oskar_output/sim.ms",
            "channel_bandwidth_hz": 5e4,
            "time_average_sec": 10.0,
            "uv_filter_max": 1000,
            "uv_filter_units": "Wavelengths",
        },
        "sky": {"oskar_sky_model/file": "BTA/sky_model.osm"},
    }

    DEFAULT_IMAGER_SETTINGS: dict = {
        "General": {"app": "oskar_imager"},
        "image": {
            "use_gpus": True,
            "channel_snapshots": True,
            "input_vis_data": "BTA/oskar_output/sim.ms",
            "root_path": "BTA/oskar_output/sim_image",
            # Below values are modified dynamically =>
            "fov_deg": 1.5,
            "size": 100,
            # <=
        },
    }

    PRIMARY_GENERAL_SETTINGS: dict = {"General": {}}

    DEFAULT_GENERAL_SETTINGS: dict = (
        DEFAULT_IMAGER_SETTINGS | DEFAULT_INTERFEROMETER_SETTINGS | PRIMARY_GENERAL_SETTINGS
    )

    # Legal settings keywords
    LEGAL_INTERFEROMETER_HEADINGS: set[str] = {
        "simulator",
        "sky",
        "telescope",
        "observation",
        "interferometer",
    }
    LEGAL_IMAGER_HEADINGS: set[str] = {"image"}

    # Load yuxiang's h5 data
    # Properties: size = (400, 400, 400) px; voxels = (1.5, 1.5, 1.5) cMPc; z_ref = ~7 (box #1), ~8 (box #2)
    @staticmethod
    def load_coeval_templates(template_switch: bool, oskar_parent_dir: str = "~") -> np.ndarray:
        """Loads the values for the coeval templates.

        :param template_switch: If true load template 1, if false load template 2.
        :param oskar_parent_dir: The directory containing the oskareor.data folder (default is the home folder).

        :return coeval_bt_array: The brightness temperature array corresponding to one of two coeval boxes.
        """

        try:
            coeval_template = h5py.File(
                ofmg.expand_path(
                    oskar_parent_dir
                    + "/oskareor.data/simulations/legacy_templates/yuxiang"
                    + ("1" if template_switch else "2")
                    + ".h5"
                ),
                "r",
            )
        except FileNotFoundError:
            print("WARNING: The coeval template h5 file was not found. Random data values will be used instead.")

            return np.random.rand(400, 400, 400)

        return np.array(coeval_template.get("BrightnessTemp")["brightness_temp"])

    @staticmethod
    def select_option(options: dict, selection: str) -> dict:
        """Select an option from an option dictionary.

        :param options (dict): The option dictionary. Dictionary must have a value structure of (option synonyms, description, *other).
        :param selection (str): The option to select. If empty return the whole option dictionary.

        :return option_dict (dict): A dictionary containing only the option or the entire option dictionary
        """

        if selection == "":
            return options

        for option in ostr.recase_iterable(set(options.keys())):
            if selection.lower() == option or selection.lower() in options[option][0]:
                return {option: options[option]}

        raise ValueError("Option " + selection + " is not a valid option.")

    @staticmethod
    def display_options(options: dict, selection: str = "", print_options: bool = True) -> dict:
        """Select and/or display option(s) from an option dictionary.

        :param options (dict): The option dictionary. Dictionary must have a value structure of (option synonyms, description, *other).
        :param selection (str): The option to select. If empty return the whole option dictionary.
        :param print_options (bool): If true print the help for the option dictionary. If false return dictionary only.

        :return options_dict (dict): A dictionary containing only the option or the entire option dictionary.
        """

        options = OSKARHelper.select_option(options, selection)

        if print_options:
            print("============")
            for option in options:
                print(
                    "OPTION: "
                    + option.upper()
                    + "\n"
                    + "Description: "
                    + options[option][1]
                    + "\n"
                    + "Synonyms: "
                    + option
                    + ", "
                    + ", ".join(options[option][0])
                    + "\n"
                    + "Preset: "
                    + option
                    + "\n"
                    + "============"
                )
            print("Option keys will ignore case.")
            print("============")

        return options


class FunctionTimer:
    """A wrapper object for creating a function with a timer."""

    def __init__(self, func, msg=None):
        self.start_time = None
        self.end_time = None
        self.func = func
        self.msg = msg

    def execute(self, *func_args: tuple, **func_kwargs: dict[str, Any]) -> Any:
        """Executes a function and times its execution, prints the elapsed time to terminal.

        :param *func_args (tuple): Positional arguments for the function.
        :param **func_kwargs (tuple): Keyword arguments for the function.

        :return output (Any): Whatever is returned from the base function
        """

        self.start_time = timer()

        output = self.func(*func_args, **func_kwargs)

        self.end_time = timer()

        if self.msg is not None:
            self.reprint(self.msg)

        return output

    def reprint(self, msg):
        """Print the time elapsed with a custom message."""
        print(msg, "\nTime elapsed:", f"{self.get_time():4.3f}", "seconds.")

    def get_time(self) -> float:
        """Get the elapsed time."""
        if self.end_time is None or self.start_time is None:
            raise RuntimeError("Function has not been called yet.")
        else:
            return self.end_time - self.start_time

    @staticmethod
    def execute_and_time(func: Callable, *func_args: tuple, msg: str = None, **func_kwargs: dict[str, Any]) -> Any:
        """Executes a function and times its execution, prints the elapsed time to terminal.

        :param func (Callable): The function in question to call.
        :param msg (str): The string message to print to terminal on its own line.
        :param *func_args (tuple): Positional arguments for the function.
        :param **func_kwargs (tuple): Keyword arguments for the function.

        :return output (Any): Whatever is returned from the base function
        """

        func_timer = FunctionTimer(func, msg)

        return func_timer.execute(*func_args, **func_kwargs)


# TODO: Create oskareor.data dir creator
