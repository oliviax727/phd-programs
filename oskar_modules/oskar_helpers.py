"""
The oskar_helpers module contains a series of helper functions and variables for the entire simulation-to-power-spectra pipeline.
"""

# System imports
import os
import configparser as cfp

# Mathematics and calculations
import astropy.constants as c
import astropy.units as u
import numpy as np
import h5py

# Data handling and statistics
from scipy.stats import circmean as circmean_radians

# Astropy extras
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time, TimeDelta
from astropy.cosmology import FlatLambdaCDM as fmodel
from astropy.cosmology import z_at_value as getz

# Reformat helper functions
class OSKARHelper():
    """
        Helper functions and constants specific to the Reformatding and oskar handling process.
    """

    # Define Constants
    SKA_REF_LOC = EarthLocation.of_site("SKA-LOW")
    OBS_LEN_4HR = TimeDelta(4 * u.hr)
    REF_TIME    = Time(val="2025-03-03T05:30:00.00", format='isot', scale='utc')
    ZENITH_530  = SkyCoord(ra=0*u.deg, dec=-27*u.deg, frame='icrs') # SKA-Low Zenith at 5:30 am 2025-03-03
    ZERO_RADEC  = SkyCoord(ra=0*u.deg, dec=0*u.deg, frame='icrs') # Centre RA/dec
    OSKAR_SIF   = "~/.oskar/OSKAR-2.12.2-Python3.sif"
    OSKAR_BIN   = "~/.oskar/bin/"
    TELESCOPE   = "~/.oskar/SKA-Low_telescope_models/SKA-Low_AAstar_original_rigid-rotation.tm"
    SIGMA_F     = (np.sqrt(c.k_B / (1.008 * c.u * (21.106 * u.cm)**2))).to(u.Hz*u.K**-0.5).value

    # Define default settings
    DEFAULT_INTERFEROMETER_SETTINGS = {
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

    DEFAULT_IMAGER_SETTINGS = {
        "General": {
            "app": "oskar_imager"
        },
        "image": {
            "use_gpus": False,
            "channel_snapshots": "false",
            "input_vis_data": "output/sim.ms",
            "root_path": "output/sim_image"
        }
    }

    # Calculated settings: [observation] start_frequency_hz, num_channels, frequency_inc_hz, phase_centre_ra_deg, phase_centre_dec_deg, length, start_time_utc
    # Calculated settings: [image] fov_deg, size
    DEFAULT_DYNAMIC_GENERAL_SETTINGS = {
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

    PRIMARY_GENERAL_SETTINGS = { "General": {} }

    DEFAULT_GENERAL_SETTINGS = DEFAULT_IMAGER_SETTINGS | DEFAULT_INTERFEROMETER_SETTINGS | DEFAULT_DYNAMIC_GENERAL_SETTINGS | PRIMARY_GENERAL_SETTINGS

    # Legal settings keywords
    LEGAL_INTERFEROMETER_HEADINGS = { "simulator", "sky", "telescope", "observation", "interferometer"}
    LEGAL_IMAGER_HEADINGS = { "image" }

    # Load yuxiang's h5 data
    # Properties: size = (400, 400, 400) px; voxels = (1.5, 1.5, 1.5) cMPc; z_ref = ~7 (box #1), ~8 (box #2)
    COEVAL_TEMPLATE_1        = h5py.File('/home/olivia/.oskar/simulations/legacy_templates/yuxiang1.h5', 'r')
    COEVAL_TEMPLATE_2        = h5py.File('/home/olivia/.oskar/simulations/legacy_templates/yuxiang2.h5', 'r')
    COEVAL_TEMPLATE_VALUES_1 = np.array(COEVAL_TEMPLATE_1.get('BrightnessTemp')['brightness_temp'])
    COEVAL_TEMPLATE_VALUES_2 = np.array(COEVAL_TEMPLATE_2.get('BrightnessTemp')['brightness_temp'])

    @staticmethod
    def select_option(options, selection):
        """
        Select an option from an option dictionary.

        :param options: The option dictionary. Dictionary must have a value structure of (option synonyms, description, *other).
        :param selection: The option to select. If empty return the whole option dictionary.

        :return: A dictionary containing only the option or the entire option dictionary
        """

        if selection == "":
            return options

        for option in options:
            if selection.lower() == option or selection.lower() in options[option][0]:
                return { option : options[option] }
            
        raise ValueError("Option "+selection+" is not a valid option.")

    @staticmethod
    def display_options(options, print_options=True, selection=""):
        """
        Select and/or display option(s) from an option dictionary.

        :param options: The option dictionary. Dictionary must have a value structure of (option synonyms, description, *other).
        :param selection: The option to select. If empty return the whole option dictionary.
        :param print_options: If true print the help for the option dictionary. If false return dictionary only.

        :return: A dictionary containing only the option or the entire option dictionary.
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

    @staticmethod
    def expand_path(path):
        """ Expands a filepath to produce an absolute path. """
        return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))

    @staticmethod
    def dir_list_sorted(dir_):
        """ Retreives a frequency-sorted list of OSM or FITS file names from a given directory. """

        files = np.array(os.listdir(dir_))

        for file in files:
            if "_osm.fits" in file:
                files = files[files != file]

        sort_type = [('file', 'O'), ('num', int)]

        sorted_files, _ = zip(*np.sort(np.array(list(map(
                lambda x: (x, int(x.split('.')[1].split("_")[0])),
                files
                )), dtype=sort_type), order="num"))

        return sorted_files
    
    @staticmethod
    def find_replace_line(file_name, find_line, replace_line):
        """
        Replace a given line in a settings.ini file given the line is equal to a special string.
        
        :param file: The file to perform the find-and-replace.
        :param find_line: The special keyword to trigger a replace.
        :param replace_line: The line to replace the preset.
        """

        with open(file_name, 'r', encoding='utf-8') as file:
            lines = file.readlines()

        with open(file_name, 'w', encoding='utf-8') as file:
            for line in lines:
                if line.startswith(find_line):
                    file.write(replace_line+"\n")
                else:
                    file.write(line)

    @staticmethod
    def read_settings_to_dictionary(file):
        """
        Convert a settings ini file into a dictionary using configparser.

        :param file: The specific ini file to read.

        :return: The settings in the form of a dictionary.
        """

        config = cfp.ConfigParser()
        config.read_file(open(file, encoding='utf-8'))

        return { s: dict(config.items(s)) for s in config.sections() }
    
    @staticmethod
    def save_settings_from_dictionary(file, settings_dict):
        """
        Save a dictionary contents as an ini file using configparser.

        :param file: The specific ini file to save the data to.
        """

        config = cfp.ConfigParser()
        config.read_dict(settings_dict)

        with open(file, 'w', encoding='utf-8') as settings_file:
            config.write(settings_file)

class Maths():
    """
    Mathematical helper functions for various calculations used throughout the document. All angle units are in radians unless if otherwise specified.
    """

    # Angular distance calculations in degrees
    @staticmethod
    def normalise_angle(t):
        """ Calculate the normalised angle """
        return (t % 360 + 180) % 360 - 180
    
    normalise = normalise_angle
    norm = normalise_angle

    @staticmethod
    def denormalise_angle(t):
        """ Denormalise an angle i.e. convert to [0, 360). """
        return (t % 360 + 360) % 360
    
    denormalise = denormalise_angle
    denorm = denormalise_angle

    @staticmethod
    def angle_delta(a, b):
        """ Calculate the unnormalised difference between two angles. """
        return Maths.norm(a)-Maths.norm(b)

    delta = angle_delta

    @staticmethod
    def angle_difference(a, b):
        """ Calculate the true difference between two angles. """
        return min(abs(Maths.delta(a, b)),360-abs(Maths.delta(a, b)))

    angle_difference = np.vectorize(angle_difference)
    diff = angle_difference

    @staticmethod
    def signed_angle_difference(a, b):
        """ Calculate the signed difference between two angles. """
        return Maths.norm(Maths.delta(a, b) - 360) if (Maths.delta(a, b)) > 180 else ((Maths.delta(a, b) + 360) if Maths.delta(a, b) < -180 else Maths.delta(a, b))

    signed_angle_difference = np.vectorize(signed_angle_difference)
    diff_sgn = signed_angle_difference

    # Scipy circmean function for degrees
    @staticmethod
    def circmean_deg(angles, bounds=(0, 360)):
        """ Calculate the mean using modular arithmetic, a modified version of scipy.stats.circmean. """
        return np.rad2deg(circmean_radians(np.deg2rad(angles), high=np.deg2rad(bounds[1]), low=np.deg2rad(bounds[0])))

    circmean = circmean_deg

    # Normal, sinusoid, and sinc functions for convenience
    @staticmethod
    def normal(x, mean=0, var=1, amp=1/np.sqrt(2*np.pi)):
        """ A simple normal distribution function. """
        return amp * np.exp(-(x-mean)**2/(2*var))/np.sqrt(2*np.pi*var)
    
    gaussian = normal
    
    @staticmethod
    def sinusoid(x, f=1, ph=0, amp=1):
        """ A simple sinusoid function. """
        return amp * np.cos(2*np.pi*f*x+ph)

    @staticmethod
    def sinc(x, f=1, ph=0, amp=1):
        """ A simple sinc function. """
        return amp * np.nan_to_num(np.sin(2*np.pi*f*x+ph)/(2*np.pi*f*x+ph), nan=1, posinf=1, neginf=1)
    
    # Convert FWHM to Variance and vice versa
    @staticmethod
    def full_width_at_half_maximum(x):
        """ Calculate the FWHM from the standard deviation. Distribution is assumed as gaussian. """
        return x * 2 * np.sqrt(2*np.log(2))
    
    FWHM = full_width_at_half_maximum

    @staticmethod
    def standard_deviation(x):
        """ Calculate the standard deviation from the FWHM. Distribution is assumed as gaussian. """
        return x / (2 * np.sqrt(2*np.log(2)))
    
    STDEV = standard_deviation

    # l, m, n to RA, dec
    @staticmethod
    def lm_to_radec(l, m, phase_centre=OSKARHelper.ZERO_RADEC):
        """ Convert the l, m plane coordinates to RA and dec. """
        d0 = phase_centre.dec.to_value(u.rad)
        a0 = phase_centre.ra.to_value(u.rad)
        n = np.sqrt(1-l**2-m**2)

        dec = np.arcsin(m*np.cos(d0)+n*np.sin(d0))
        ra = a0 + np.arctan(l/(n*np.cos(d0)-m*np.sin(d0)))

        return [ra, dec]

class Cosmo():
    """
    Defines a specific cosmological model for other classes to refer to. Assumes a flat universe.

    :param h0: The Hubble constant.
    :param omega_m_0: The dimensionless matter density.
    :param omega_b_0: The dimensionless baryonic matter density.
    :param cosmo: The cosmological ΛCDM model.
    """

    def __init__(self, h0=100, omega_m_0=0.31, omega_b_0=0.048):
        self.h0        = h0 * u.km / u.s / u.Mpc # Set Hubble Constant to 100 h, with h being dimensionless hubble parameter
        self.omega_m_0 = omega_m_0
        self.omega_b_0 = omega_b_0
        self.cosmo     = fmodel(h0=h0, omega_m_0=omega_m_0, omega_b_0=omega_b_0) # Flat ΛCDM means Dark Energy density is 0.69

    # Redshift to comoving distance
    def z_to_dz(self, z):
        """ Convert redshift to comoving distance. """
        return self.cosmo.comoving_distance(z)

    # Comoving distance to redshift
    def dz_to_z(self, dz):
        """ Convert comoving distance to redshift. """
        return getz(self.cosmo.comoving_distance, dz)

    # Redshift to frequency in GHz
    @staticmethod
    def z_to_f(z, f0 = 1.420e9 * u.Hz):
        """ Convert redshift to frequency. By default it assumes that the start frequency is the 21 cm line. """
        return f0 / (z + 1)
    
    # Frequency in GHz to redshift
    @staticmethod
    def f_to_z(f, f0 = 1.420e9 * u.Hz):
        """ Convert frequency to redshift. By default it assumes that the start frequency is the 21 cm line. """
        return (f0 / f) - 1
