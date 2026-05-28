"""
The oskar_helpers module contains a series of helper functions and variables for simulation-to-power-spectra pipeline. Requires an existing and robust .oskar file in the home directory.
"""

# Mathematics and calculations
import h5py
import numpy as np

from oskareor.skalow_calc import OSKARFileConfig as ofc

# OSKAR EOR helper class - dependant on ~/.oskar
class OSKARHelper():
    """
    Helper functions and constants specific to handling OSKAR and the Reformatter.
    """

    # Define Constants
    OSKAR_SIF: str = "~/.oskar/OSKAR-2.12.2-Python3.sif"
    OSKAR_BIN: str = "~/.oskar/bin/"
    TELESCOPE: str = "~/.oskar/SKA-Low_telescope_models/SKA-Low_AAstar_original_rigid-rotation.tm"

    # Define default settings
    DEFAULT_INTERFEROMETER_SETTINGS: dict = {
        "General": {
            "app": "oskar_sim_interferometer"
        },
        "simulator": {
            "use_gpus": False
        },
        "observation" : {
            "num_time_steps": 24,
        },
        "telescope": {
            "input_directory": "telescope_model",
            "apeture_array/element_pattern/enable_numerical": False,
            "apeture_array/array_pattern/element/x_gain": 1.0,
            "apeture_array/array_pattern/element/y_gain": 1.0,
            "apeture_array/array_pattern/element/x_gain_error_time": 0.0015057,
            "apeture_array/array_pattern/element/y_gain_error_time": 0.0015057,
            "apeture_array/array_pattern/element/x_phase_error_fixed_deg": 0.0,
            "apeture_array/array_pattern/element/y_phase_error_fixed_deg": 0.0,
            "apeture_array/array_pattern/element/x_phase_error_time_deg": 0.0015057,
            "apeture_array/array_pattern/element/y_phase_error_time_deg": 0.0015057
        },
        "interferometer": {
            "oskar_vis_filename": "oskar_output/vis.vis",
            "ms_filename": "oskar_output/sim.ms",
            "channel_bandwidth_hz": 5e4,
            "time_average_sec": 10.0,
            "uv_filter_max": 1000,
            "uv_filter_units": "Wavelengths"
        },
        "sky": {}
    }

    DEFAULT_IMAGER_SETTINGS: dict = {
        "General": {
            "app": "oskar_imager"
        },
        "image": {
            "use_gpus": False,
            "channel_snapshots": "false",
            "input_vis_data": "oskar_output/sim.ms",
            "root_path": "oskar_output/sim_image"
        }
    }

    # Calculated settings: [observation] start_frequency_hz, num_channels, frequency_inc_hz, phase_centre_ra_deg, phase_centre_dec_deg, length, start_time_utc
    # Calculated settings: [image] fov_deg, size
    DEFAULT_DYNAMIC_GENERAL_SETTINGS: dict = {
        "observation" : {
            "start_frequency_hz": 200e6,
            "num_channels": 100,
            "frequency_inc_hz": 140e3,
            "phase_centre_ra_deg": 0.0,
            "phase_centre_dec_deg": -27.0,
            "length": "4:00:00.00",
            "start_time_utc": "2025-03-03 03:30:00.00"
        },
        "image" : {
            "fov_deg": 1.5,
            "size": 100
        }
    }

    PRIMARY_GENERAL_SETTINGS: dict = { "General": {} }

    DEFAULT_GENERAL_SETTINGS: dict = DEFAULT_IMAGER_SETTINGS | DEFAULT_INTERFEROMETER_SETTINGS | DEFAULT_DYNAMIC_GENERAL_SETTINGS | PRIMARY_GENERAL_SETTINGS

    # Legal settings keywords
    LEGAL_INTERFEROMETER_HEADINGS: set   = { "simulator", "sky", "telescope", "observation", "interferometer"}
    LEGAL_IMAGER_HEADINGS: set           = { "image" }

    # Load yuxiang's h5 data
    # Properties: size = (400, 400, 400) px; voxels = (1.5, 1.5, 1.5) cMPc; z_ref = ~7 (box #1), ~8 (box #2)
    COEVAL_TEMPLATE_1: h5py.File            = h5py.File(ofc.expand_path('~/.oskar/simulations/legacy_templates/yuxiang1.h5'), 'r')
    COEVAL_TEMPLATE_2: h5py.File            = h5py.File(ofc.expand_path('~/.oskar/simulations/legacy_templates/yuxiang2.h5'), 'r')
    COEVAL_TEMPLATE_VALUES_1: np.ndarray    = np.array(COEVAL_TEMPLATE_1.get('BrightnessTemp')['brightness_temp'])
    COEVAL_TEMPLATE_VALUES_2: np.ndarray    = np.array(COEVAL_TEMPLATE_2.get('BrightnessTemp')['brightness_temp'])

    @staticmethod
    def select_option(options: dict, selection: str) -> dict:
        """
        Select an option from an option dictionary.

        :param options (dict): The option dictionary. Dictionary must have a value structure of (option synonyms, description, *other).
        :param selection (str): The option to select. If empty return the whole option dictionary.

        :return option_dict (dict): A dictionary containing only the option or the entire option dictionary
        """

        if selection == "":
            return options

        for option in ofc.recase_iterable(set(options.keys())):
            if selection.lower() == option or selection.lower() in options[option][0]:
                return { option : options[option] }
          
        raise ValueError("Option "+selection+" is not a valid option.")

    @staticmethod
    def display_options(options: dict, selection: str = "", print_options: bool = True) -> dict:
        """
        Select and/or display option(s) from an option dictionary.

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
                    "OPTION: "      + option.upper()                                + "\n" +
                    "Description: " + options[option][1]                            + "\n" +
                    "Synonyms: "    + option + ", " + ", ".join(options[option][0]) + "\n" +
                    "Preset: "      + option                                        + "\n" +
                    "============"
                    )
            print("Option keys will ignore case.")
            print("============")
          
        return options
