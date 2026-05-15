#!.venv/bin/python
"""
The oskar_regrid module contains all funtions that help convert early universe simulations into sky models, measurement sets, and power spectra for the analysis of simulated SKA observations.
"""

# Import the stuffs
import os
import warnings
import subprocess
import configparser as cfp
import astropy.constants as c
import astropy.units as u
import h5py
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time, TimeDelta
from astropy.cosmology import FlatLambdaCDM as fmodel
from astropy.cosmology import z_at_value as getz
from scipy.interpolate import make_interp_spline as misp
from scipy.stats import circmean as circmean_radians

# TODO: Turn into pip project (later)

# TODO: Clean pylint errors
# pylint: disable=invalid-name

# pylint: disable-next=unused-import
from matplotlib import pyplot as plt

# Regrid helper functions
class RegridHelper():
    """
        Helper functions and constants specific to the regridding and oskar handling process.
    """

    # Define Constants
    SKA_REF_LOC = EarthLocation.of_site("SKA-LOW")
    OBS_LEN_4HR = TimeDelta(4 * u.hr)
    REF_TIME    = Time(val="2025-03-03T05:30:00.00", format='isot', scale='utc')
    ZENITH_530  = SkyCoord(ra=0*u.deg, dec=-27*u.deg, frame='icrs') # SKA-Low Zenith at 5:30 am 2025-03-03
    ZERO_RADEC  = SkyCoord(ra=0*u.deg, dec=0*u.deg, frame='icrs') # Centre RA/Dec
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

        options = RegridHelper.select_option(options, selection)

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
        return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))

    @staticmethod
    def dir_list_sorted(dir_):
        """
        Retreives a frequency-sorted list of OSM or FITS file names from a given directory.

        :param dir_: The directory to pull the files from.
        :return: The sorted list of files.
        """

        files = np.array(os.listdir(dir_))

        for file in files:
            if "_osm.fits" in file:
                files = files[files != file]

        sort_type = [('file', 'O'), ('num', int)]
        sort_prep = lambda x: (x, int(x.split('.')[1].split("_")[0]))

        sorted_files, _ = zip(*np.sort(np.array(list(map(sort_prep, files)), dtype=sort_type), order="num"))

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

    # l, m, n to RA, Dec
    @staticmethod
    def lm_to_radec(l, m, phase_centre=RegridHelper.ZERO_RADEC):
        d0 = phase_centre.dec.to_value(u.rad)
        a0 = phase_centre.ra.to_value(u.rad)
        n = np.sqrt(1-l**2-m**2)

        Dec = np.arcsin(m*np.cos(d0)+n*np.sin(d0))
        Ra = a0 + np.arctan(l/(n*np.cos(d0)-m*np.sin(d0)))

        return [Ra, Dec]

class Cosmo():
    """
    Define a static cosmology method for other classes to use.
    """

    def __init__(self, H0=100, Om0=0.31, Ob0=0.048):
        self.H0    = H0 * u.km / u.s / u.Mpc # Set Hubble Constant to 100 h, with h being dimensionless hubble parameter
        self.Om0   = Om0
        self.Ob0   = Ob0
        self.cosmo = fmodel(H0=H0, Om0=Om0, Ob0=Ob0) # Flat ΛCDM means Dark Energy density is 0.69

    # Redshift to comoving distance
    def z_to_Dz(self, z): return self.cosmo.comoving_distance(z)

    # Comoving distance to redshift
    def Dz_to_z(self, Dz): return getz(self.cosmo.comoving_distance, Dz)

    # Redshift to frequency in GHz
    @staticmethod
    def z_to_f(z): return 1.42e9 * u.Hz / (z + 1)

class Regrid():
    """
    The regridding class contains functions relating to translating simulation data to OSKAR output data.
    """

    @staticmethod
    def brightness_temperature_to_flux(Tb, fxy, dtheta, dphi):
        # Calculate pixelated luminosity
        Fv = 2 * c.k_B.value * fxy**2 * Tb * (dtheta * dphi) / c.c.value ** 2
        Fv = Fv * 1e26 # Converts to Jansky
        return Fv
    
    @staticmethod
    def brightness_temperature_to_linewidth(Tb):
        # Calculate linewidth in Hz
        return (Tb ** 0.5) * RegridHelper.SIGMA_F
    
    TEMPLATE_PRESETS = {
        "gaussian" : (
            { "normal", "gauss", "bell", "n", "g", "b" },
            "Rotationally symmetric centered gaussian plane with FWHM = d(t)/3.",
            lambda p: Maths.gaussian(p['r'], var=Maths.STDEV(p['d'][2]/3)**2, amp=p['T_max'])
            ),
        "flat"     : (
            { "plane", "constant", "const", "c", "f" },
            "Constant temperature for every pixel.",
            lambda p: p['T_max']
            ),
        "random"   : (
            { "rand", "r" },
            "Random values for every cell from 0 to the defined scale (T_max).",
            lambda p: p['T_max'] * np.random.rand()
            ),
        "sinusoid" : (
            { "sinc", "interference", "fringe", "i", "intf", "s" },
            "Rotationally symmetric centered sinc function with freq = 1/d(t).",
            lambda p: np.abs(Maths.sinc(p['r'], f=1/p['d'][2], amp=p['T_max']))
            ),
        "point"    : (
            { "delta", "source", "p" },
            "A point source in the direct centre of the field.",
            lambda p: p['T_max'] if (p['i'] == p['d'][0] // 2 and p['j'] == p['d'][1] // 2) else 0
            ),
        "dark"     : (
            { "clear", "empty", "d" },
            "A completely clear sky.",
            lambda p: 0,
            ),
        "coeval1"  : (
            { "coeval 1", "1", "yuxiang1", "yuxiang 1", "y1", "c1" },
            "One of two simulation boxes, cocentric with the desired values box, the original model has d = (400, 400, 400).\n"
            + "If d(a) < 400 then the box outer edges will be cropped and if d(a) > 400 the box will repeat beyond 400 px from the centre.",
            lambda p: RegridHelper.COEVAL_TEMPLATE_VALUES_1[*((200 + ((np.array([p['i'], p['j'], p['t']]) - (np.array(p['d']) // 2)))) % 400)]
            ),
        "coeval2"  : (
            { "coeval 2", "2", "yuxiang2", "yuxiang 2", "y2", "c2" },
            "The second of two simulation boxes, cocentric with the desired values box, the original model has d = (400, 400, 400).\n"
            + "If d(a) < 400 then the box outer edges will be cropped and if d(a) > 400 the box will repeat beyond 400 px from the centre.",
            lambda p: RegridHelper.COEVAL_TEMPLATE_VALUES_2[*((200 + ((np.array([p['i'], p['j'], p['t']]) - (np.array(p['d']) // 2)))) % 400)]
            )
    }

    @staticmethod
    def display_template_presets(print_presets=True, filter_preset=""):
        """
        Return and/or display all available templates, their names, and their descriptions. All templates except random, coeval1, and coeval2 are identical in all t-dimension voxels.

        :param print_presets: If true print the templates to console and return the dictionary. If false only return the dictionary.
        :param filter_preset: If true print/return the information for only one specific template entry, according to the given string.
        
        :return: The dictionary containing all available templates or a dictionary of the specific desired template.
        """

        return RegridHelper.display_options(Regrid.TEMPLATE_PRESETS, print_options=print_presets, selection=filter_preset)

    @staticmethod
    def mock_values(preset, scale = 10, d = (100, 100, 100), special = None):
        """
        Create an array of mock simulation values.

        :param preset: Mock brightness temperature array format. Run Regrid.display_template_presets for more information.
        :param scale: a.k.a. `T_max`. The maximum Kelvin value for the whole array, acts as a normalisation factor.
        :param d: The size of the values datacube.
        :param special: A custom lambda function that takes the dictionary of parameters (`d`, `i`, `j`, `x`, `y`, `t`, `r`, `T_max`) and returns a float, treat the preset parameter as a custom name.
        
        Note that `x` and `y` are positioned so that the centermost pixel is (0, 0) whereas `i` and `j` are the standard array values array indicies. Only `d`, `i`, and `j` are indicies, the others should be treated as floats.

        :return: Mock brightness temperature values.
        """

        # Select specific template
        if special is None:
            selection = Regrid.display_template_presets(False, preset)
            preset = list(selection.keys())[0]
            func = selection[preset][2]
        else:
            func = special

        # Create initial array
        print("Creating mock brightness temperatures for the template: "+preset)
        values = np.zeros(d).astype(np.float64)

        # Iterate through all elements of the dictionary
        for t in range(d[2]):
            for i in range(d[0]):
                for j in range(d[1]):
                    # Define the iteration dictionary
                    params = {
                        "d" : d, "i" : i, "j" : j, "t": t,
                        "x" : i - d[0]/2, "y" : j - d[1]/2,
                        "r": np.sqrt((i - d[0]/2)**2 + (j - d[1]/2)**2),
                        "T_max" : scale
                        }
                    
                    # Populate array cell
                    values[i, j, t] = func(params)

            print("\rTime step #", t, end="")
        
        print("\nValues created!")

        return values.astype(np.float64)

    @staticmethod
    def convert_H5_coeval_to_csv(h5_location, save_data=False, outdir='', name="out_h5_data"):
        """
        Extract h5 data from Yuxiang Qin's simulations and either return the data objects or save to a seperate CSV.

        :param h5_location: The file location of the h5 data.
        :param save_data: If true, output data to a CSV and text file, specified by the outdir parameter.
        :param outdir: The directory to output both CSV and text information.
        :param name: The file name template to be saved to.

        :return: The numpy values array in Kelvin, the shape of the array, the refrence redshift, the voxel size in Mpc, and the simulation box cosmology.
        """

        file = h5py.File(h5_location, 'r')

        # Get BT data
        bt_data = np.array(file.get('BrightnessTemp')['brightness_temp'])

        # Get Box and Voxel dimensions
        box_len = file.get('user_params').attrs['BOX_LEN']#* file.get('cosmo_params').attrs['hlittle']
        vox = np.ones(3) * box_len / bt_data.shape[0]

        # Define cosmology with H0=100h
        cosmology = Cosmo(
            Om0 = file.get('cosmo_params').attrs['OMm'],
            Ob0 = file.get('cosmo_params').attrs['OMb']
        )

        # Transform intitial redshift
        z_mid = file.attrs['redshift']
        z_ref = cosmology.Dz_to_z(cosmology.z_to_Dz(z_mid) - u.Mpc*box_len/2)

        if save_data:
            np.savetxt(outdir+'/'+name+'.csv', bt_data, delimiter=", ")
            np.savetxt(outdir+'/'+name+'.txt', np.array([bt_data.shape, z_ref, vox, cosmology]), delimiter=", ")

        return bt_data, bt_data.shape, z_ref, vox, cosmology
    
    @staticmethod
    def convert_H5_lightcone_to_csv(h5_location, save_data=False, outdir='', name="out_h5_data"):
        """
        Extract h5 data from Yuxiang Qin's lightcone simulations and either return the data objects or save to a seperate CSV.

        :param h5_location: The file location of the h5 data.
        :param save_data: If true, output data to a CSV and text file, specified by the outdir parameter.
        :param outdir: The directory to output both CSV and text information.
        :param name: The file name template to be saved to.
        :return: The numpy values array in Kelvin, the shape of the array, the refrence redshift, the voxel size in Mpc, and the simulation box cosmology.
        """

        file = h5py.File(h5_location, 'r')

        # Get BT data
        bt_data = np.array(file.get('lightcones/brightness_temp'))

        # Define cosmology with H0=100h
        cosmology = Cosmo(
            Om0 = file.get('cosmo_params').attrs['OMm'],
            Ob0 = file.get('cosmo_params').attrs['OMb']
        )

        # Get Box and Voxel dimensions
        box_len = file.get('user_params').attrs['BOX_LEN'] * file.get('cosmo_params').attrs['hlittle']
        rs = np.array(list(file.get('node_redshifts')))
        vox = np.array([box_len, box_len, abs(cosmology.z_to_Dz(rs[0])-cosmology.z_to_Dz(rs[-1])).to_value(u.Mpc)]) / bt_data.shape
        
        # Transform intitial redshift
        z_ref = np.min(rs)

        if save_data:
            np.savetxt(outdir+'/'+name+'.csv', bt_data, delimiter=", ")
            np.savetxt(outdir+'/'+name+'.txt', np.array([bt_data.shape, z_ref, vox, cosmology]), delimiter=", ")

        return bt_data, bt_data.shape, z_ref, vox, cosmology
    
    @staticmethod
    def convert_H5_to_csv(h5_location, save_data=False, outdir='', name="out_h5_data", coeval=True):
        """
        Extract h5 data from Yuxiang Qin's simulations and either return the data objects or save to a seperate CSV.

        :param h5_location: The file location of the h5 data.
        :param save_data: If true, output data to a CSV and text file, specified by the outdir parameter.
        :param outdir: The directory to output both CSV and text information.
        :param name: The file name template to be saved to.
        :param coeval: whether or not the box is coeval or lightcone.
        :return: The numpy values array in Kelvin, the shape of the array, the refrence redshift, the voxel size in Mpc, and the simulation box cosmology.
        """

        if coeval:
            return Regrid.convert_H5_coeval_to_csv(h5_location=h5_location, save_data=save_data, outdir=outdir, name=name)
        else:
            return Regrid.convert_H5_lightcone_to_csv(h5_location=h5_location, save_data=save_data, outdir=outdir, name=name)
        
    @staticmethod
    def transform_datacube_units(values, voxels, z_ref = 7, require_regrid = True, max_freq_res = 100 * u.MHz, cosmology=Cosmo()):
        """
        Transform a datacube with dimensions x, y, t (cMpc x cMpc x cMpc) to ⍺, δ, f (rad x rad x Hz),

        :param values: The simulation datacube in units of Kelvin.
        :param voxels: An array describing a series of voxel dimensions corresponding to each simulation datacube voxel element in units of (cMpc, cMpc, cMps).
        :param z_ref: Refrence redshift, the ending redshift of the simulation.
        :param require_regrid: If true then always regrid frequency bins, if false, regrid only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution.
        :param v: If all voxels are the same, provides the initial voxel dimensions in cMpc in dimensions (x, y, t), and auto-generates the voxel configuration array.
        :param osm_output: The relative path to save the osm file to.
        :param cosmology: The specific cosmology parameters in the form of a custom Cosmo object.
        """

        # Configure d variable
        d = np.shape(values)

        # Set regrid flag
        regrid_flag = require_regrid

        # Calculate refrence comoving dist.
        Dz_ref = cosmology.z_to_Dz(z_ref)
        f_ref = cosmology.z_to_f(z_ref)

        # Set cosmological redshift parameters
        Dz = Dz_ref
        z_prev = z_ref
        fq = f_ref.to_value(u.Hz)
        Dz_val = Dz_ref.to_value(u.Mpc)

        # Set maximum frequency resolution
        max_freq_res_hz = max_freq_res.to_value(u.Hz)

        # Set Linewidth array
        sigma_f = Regrid.brightness_temperature_to_linewidth(values)

        print("Transforming coordinates ...")
        # Main loop of creation
        for t in range(d[2]):
            for x in range(d[0]):
                for y in range(d[1]):

                    # Retreive voxel values
                    dx = voxels[x, y, t, 0]
                    dy = voxels[x, y, t, 1]
                    dt = voxels[x, y, t, 2] 
                    Tb = values[x, y, t]

                    # STEP 1 - Convert transverse comoving distances to flat angular resolution
                    Dz_pix = Dz_val+dt/2 # Alter it by the CENTRAL pixel value

                    # Calculate transformed dimensions
                    dtheta = np.arctan(dx/Dz_pix)
                    dphi = np.arctan(dy/Dz_pix)

                    # STEPS 2 & 3 - Convert line-of-sight comoving distance to frequency
                    df = 0

                    if x == y == 0:
                        # Determine corresponding redshifts
                        z_bot = z_prev
                        z_top = cosmology.Dz_to_z(Dz+dt*u.Mpc)

                        # Convert to frequency
                        f_bot = cosmology.z_to_f(z_bot)
                        f_top = cosmology.z_to_f(z_top)

                        # Store altered frequency bandwidth
                        df = np.abs((f_bot-f_top).to_value(u.Hz))

                        # STEP 7.1 - Check if regridding is needed
                        if df > max_freq_res_hz:
                            regrid_flag = regrid_flag or True

                    else:
                        df = voxels[0, 0, t, 2]

                    # STEP 4 - Convert brightness temperature to pixel flux
                    fxy = fq + df/2
                      
                    # Calculate pixelated luminosity
                    Fv = Regrid.brightness_temperature_to_flux(Tb=Tb, fxy=fxy, dtheta=dtheta, dphi=dphi)
                
                    # Save values
                    voxels[x, y, t, 0] = dtheta
                    voxels[x, y, t, 1] = dphi
                    voxels[x, y, t, 2] = df
                    values[x, y, t] = Fv
            
            # Increment distance and frequency parameters
            Dz = Dz + dt*u.Mpc # Increment the value of Dz by voxel dimension
            z_prev = z_top # Top z in current box = bottom z in next box
            Dz_val = Dz_val + dt # Increment the value of Dz by voxel dimension
            fq = fq + df # Increment cumulative frequency

            print("\rTime step #", t, end="")

        print("\nTransforming complete.")

        return values, voxels, sigma_f, f_ref, regrid_flag
        
    @staticmethod
    def regrid_datacube(values, voxels, d = None, sigma_f = None, max_freq_res=100 * u.MHz):
        """
        Regrids each spaxel of a sky model given a maximum frequency resolution.

        :param values: The simulation datacube.
        :param voxels: An array describing a series of voxel dimensions corresponding to each simulation datacube voxel element in units of (rad, rad, Hz).
        :param sigma_f: An array of same dimensions as values but containing information about the linewidth of the frequency emission profile.
        :param max_freq_res: Maximum allowable voxel frequency resolution.
        """

        print("Performing regrid ...")

        # Set maximum frequency resolution
        max_freq_res_hz = max_freq_res.to_value(u.Hz)

        # Configure d variable
        d = np.shape(values)

        for x in range(d[0]):
            for y in range(d[1]):

                print("\rSpaxel # (", x, ",", y, ")", end="")

                # Set the values as being in the middle of each bin
                freq_values = np.cumsum(voxels[x, y, :, 2]) - voxels[x, y, :, 2]/2

                # Create interpolation B-spline
                bspline = misp(freq_values, values[x, y, :])
                dspline = misp(freq_values, sigma_f[x, y, :])

                # Check maximum frequency resolution
                if np.abs(freq_values[0] - freq_values[-1])/d[2] > max_freq_res_hz:
                    # Resize d[2] to match maximum resolution
                    d[2] = np.abs(freq_values[0] - freq_values[-1])/max_freq_res_hz

                # Generate evenly-distributed frequency array
                new_freq, freq_bandw = np.linspace(freq_values[0], freq_values[-1], d[2], retstep=True)

                # Perform regridding
                new_flux = np.clip(bspline(new_freq), 0, None)
                new_sigf = np.clip(dspline(new_freq), 0, None)

                # Create array of uniform bin sizes
                freq_bins = np.ones(d[2]) * freq_bandw

                # Save variables
                voxels[x, y, :, 2] = freq_bins
                values[x, y, :] = new_flux
                sigma_f[x, y, :] = new_sigf

        print("\nRegrid complete.")
        
        return values, voxels, sigma_f
    
    @staticmethod
    def calculate_cumulative_voxels(voxels, f_ref = 200 * u.MHz, phase_ref_point = RegridHelper.ZENITH_530):
        """
        Calculate the cumulative voxel sum and centre with a refrence point and frequency.

        :param voxels: An array describing a series of voxel dimensions corresponding to each sky model datacube voxel element in units of (rad, rad, Hz).
        :param f_ref: Refrence frequency, the ending frequency of the model.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.

        :return: An array determining the specific central value of each voxel in its corresponding values array in units of (deg, deg, Hz).
        """

        # RA, Dec centering function and freq shift function
        # Add half the total and half the spaxel widths to centre the main point
        centering = (lambda x: (x - np.max(x, axis=(0, 1), keepdims=True)/2 - np.min(x, axis=(0, 1), keepdims=True)/2))

        # Add only half the cell bandwidth
        shifting = (lambda x: (x - np.min(x, axis=2, keepdims=True)/2))

        # Centre the cumulative RA and Dec sums so that the zero value is in the centre
        rasum = centering(np.cumsum(voxels[:,:,:,0], axis=0))
        decsum = centering(np.cumsum(voxels[:,:,:,1], axis=1))

        # Centre the frequency
        freqsum = f_ref.to_value(u.Hz) - shifting(np.cumsum(voxels[:,:,:,2], axis=2))

        # Calculate phase centre offsets
        source_pos = phase_ref_point.spherical_offsets_by(rasum * u.rad, decsum * u.rad)
        RAs = source_pos.ra.to_value(u.deg)
        Dcs = source_pos.dec.to_value(u.deg)

        return (RAs, Dcs, freqsum)

        
    @staticmethod
    def save_datacube_to_osm(values, voxels = None, cumulative_voxels = None, sigma_f = None, f_ref = 200 * u.MHz, phase_ref_point = RegridHelper.ZENITH_530, osm_output="regrid/osm_output/osm_output.osm"):
        """
        Saves a given datacube of flux values and voxel dimensions (RA, Dec, Freq.) to a master OSM file.

        :param values: The sky model datacube.
        :param voxels: An array describing a series of voxel dimensions corresponding to each sky model datacube voxel element (rad, rad, Hz).
        :param cumulative_voxels: A jagged array consisting of the cumulative summation of components from the voxel array (deg, deg, Hz).
        :param sigma_f: An array of same dimensions as values but containing information about the linewidth of the frequency emission profile.
        :param f_ref: Refrence frequency, the ending frequency of the model.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always regrid frequency bins, if false, regrid only when max frequency resolution is met.
        :param osm_output: The relative path to save the osm file to.
        """

        print("Configuring datacube for OSKAR file format ...")

        # Configure d variable
        d = np.shape(values)

        # Cumulative sums are more important than voxel bins now
        (RAs, Dcs, freqsum) = (None, None, None) # Keep Pylint Happy
        if cumulative_voxels is None and voxels is None:
            raise ValueError("Either an array of voxels or cumulative voxes must be provided!")
        elif voxels is None:
            (RAs, Dcs, freqsum) = cumulative_voxels
        elif cumulative_voxels is None:
            (RAs, Dcs, freqsum) = Regrid.calculate_cumulative_voxels(voxels=voxels, f_ref=f_ref, phase_ref_point=phase_ref_point)
            
        # Record data to file
        print("Recording data to .osm file")

        with open(osm_output, 'w', encoding='utf-8') as osm:
            # Clear file contents
            osm.truncate(0)

            # Add header lines
            osm.write("Format = RaD DecD I Q U V ReferenceFrequency LineWidth\n")
            osm.write("# Entries Key:\n")
            osm.write("#00.000000 +00.000000 0.0000+e00 0.0 0.0 0.0 000.000e6 0.0000+e00\n")
            osm.write("# RA       Dec        Stokes I   Q   U   V   Freq0     Linewidth\n")

            # Write OSM lines
            for x in range(d[0]):
                for y in range(d[1]):
                    for t in range(d[2]):

                        # Format data
                        rascn = np.char.zfill(np.format_float_positional(RAs[x, y, t], 6, False), 10)
                        decln = np.char.zfill(np.format_float_positional(np.abs(Dcs[x, y, t]), 6, False), 9)
                        value = np.format_float_scientific(values[x, y, t], 4, False)
                        freq0 = np.format_float_positional(freqsum[x, y, t] / 1e6, 3, False)
                        linew = np.format_float_scientific(sigma_f[x, y, t], 4, False)

                        # Add +/- value to Declinations
                        if Dcs[x, y, t] >= 0: decln = "+" + str(decln)
                        else:                 decln = "-" + str(decln)

                        # Write to OSM
                        osm.write(
                            str(rascn)    + " " + # Right Ascension
                            decln         + " " + # Declination
                            str(value)    + " " + # Intensity
                            "0.0 0.0 0.0" + " " + # Redundant Stokes Parameters
                            str(freq0)  + "e6 " + # Point source frequency
                            str(linew)    + " " + # Spectral profile linewidth
                            "\n"
                        )

                    print("\rSpaxel # (", x, ",", y, ")", end="")

        print("\nProcess complete, data saved to "+osm_output)

    # TODO: Automatically find ideal UTC time of observation
    # pylint: disable=unused-argument
    @staticmethod
    def calculate_observation_time_from_date(phase_ref_point = RegridHelper.ZENITH_530, ref_time = RegridHelper.REF_TIME, ref_location = RegridHelper.SKA_REF_LOC, observation_length = RegridHelper.OBS_LEN_4HR):
        """
        Calculates the closest ideal observation time from a given UTC date and telescope lattitude.

        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param ref_time: An astropy.time.Time object stating the desired mid-observation time.
        :param ref_location: An astropy.coordinates.EarthLocation object stating the location of the telescope on Earth.
        :param observation_length: An astropy.time.TimeDelta object that gives the length of the observation.

        :returns: The ideal UTC refrence date for the observation, and if the object goes below the horizon due to the length of observation.
        """

        obs_length_flag = False

        if obs_length_flag:
            warnings.warn("The provided observation time is too long! The object will not be in the sky for the entire duration of time.")

        return ref_time, obs_length_flag
    # pylint: enable=unused-argument

    @staticmethod
    def generate_dynamic_settings(values, voxels = None, cumulative_voxels = None, phase_ref_point = RegridHelper.ZENITH_530, f_ref = 200 * u.MHz, ref_time = RegridHelper.REF_TIME, ref_location = RegridHelper.SKA_REF_LOC, observation_length = RegridHelper.OBS_LEN_4HR, save_dynamic_settings = ""):
        """
        Generates a set of dynamically-set ini settings for OSKAR to utilise.

        :param values: The simulation datacube.
        :param voxels: An array describing a series of voxel dimensions corresponding to each sky model datacube voxel element (rad, rad, Hz).
        :param cumulative_voxels: A jagged array consisting of the cumulative summation of components from the voxel array (deg, deg, Hz).
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param f_ref: Refrence frequency, the ending frequency of the model.
        :param ref_time: An astropy.time.Time object stating the desired mid-observation time.
        :param ref_location: An astropy.coordinates.EarthLocation object stating the location of the telescope on Earth.
        :param observation_length: An astropy.time.TimeDelta object that gives the length of the observation.
        :param save_dynamic_settings: If non-empty, save the dynamic settings to an .ini file given by the path entered.

        :return: The dynamically defined settings dictionary.
        """
        # Calculated settings: [observation] start_frequency_hz, num_channels, frequency_inc_hz, phase_centre_ra_deg, phase_centre_dec_deg, length, start_time_utc
        # Calculated settings: [image] fov_deg, size

        # SETUP
        # Cumulative sums are more important than voxel bins now
        (RAs, Dcs, freqsum) = (None, None, None) # Keep Pylint Happy
        if cumulative_voxels is None and voxels is None:
            raise ValueError("Either an array of voxels or cumulative voxes must be provided!")
        elif voxels is None:
            (RAs, Dcs, freqsum) = cumulative_voxels
        elif cumulative_voxels is None:
            (RAs, Dcs, freqsum) = Regrid.calculate_cumulative_voxels(voxels=voxels, f_ref=f_ref, phase_ref_point=phase_ref_point)

        # Configure d variable
        d = np.shape(values)

        # Create deep copy of union/logical or settings set
        dynamic_settings = dict(RegridHelper.DEFAULT_GENERAL_SETTINGS)

        # BASIC CONFIG
        # Set starting frequency NB: The last channel has the lowest frequency!
        dynamic_settings['observation']['start_frequency_hz'] = np.mean(freqsum[:,:,-1])

        # Set number of channels and image size
        dynamic_settings['observation']['num_channels'] = d[2]
        dynamic_settings['image']['size'] = max(d[0], d[1])

        # Set the frequency increment
        dynamic_settings['observation']['frequency_inc_hz'] = np.mean(voxels[:,:,:,2])

        # Set phase centre RA and Dec
        dynamic_settings['observation']['phase_centre_ra_deg'] = phase_ref_point.ra.deg
        dynamic_settings['observation']['phase_centre_dec_deg'] = phase_ref_point.dec.deg

        # COORDINATE CONFIG
        # Set image field of view and size
        ref_time, _ = Regrid.calculate_observation_time_from_date(phase_ref_point=phase_ref_point, ref_time=ref_time, ref_location=ref_location, observation_length=observation_length)

        # Calculate RA dimension - use circular mean for angular averages
        rac = Maths.diff(Maths.circmean(RAs[-1,:,:]), Maths.circmean(RAs[0,:,:]))

        # Calculate Dec dimension - use circular mean for angular averages
        decc = Maths.diff(Maths.circmean(Dcs[-1,:,:], (-180, 180)), Maths.circmean(Dcs[0,:,:], (-180, 180)))

        # Set the field of view
        dynamic_settings['image']['fov_deg'] = max(rac, decc)

        # Set the observation time and length
        dynamic_settings['observation']['start_time_utc'] = (ref_time - observation_length / 2).utc.value
        dynamic_settings['observation']['length'] = str(observation_length.to_value(format='datetime'))

        # Save the dynamic settings
        if save_dynamic_settings != "":
            settings_path = RegridHelper.expand_path(save_dynamic_settings)

            RegridHelper.save_settings_from_dictionary(save_dynamic_settings, dynamic_settings)

            print("Saved dynamic and default settings to ini file: "+settings_path)
        
        return dynamic_settings

    @staticmethod
    def generate_osm_from_simulation(values, voxels = None, z_ref = 7, phase_ref_point = RegridHelper.ZENITH_530, require_regrid = True, max_freq_res = 100 * u.MHz, v = (1.5, 1.5, 1.5), osm_output="regrid/osm_output/osm_output.osm", cosmology=Cosmo(), save_dynamic_settings = "", ref_time = RegridHelper.REF_TIME, ref_location = RegridHelper.SKA_REF_LOC, observation_length = RegridHelper.OBS_LEN_4HR):
        """
        Generate a set of .osm files for an OSKAR sky model based on a Mpc**3 simulation output.

        :param values: The simulation datacube.
        :param z_ref: Refrence redshift, the ending redshift of the simulation.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always regrid frequency bins, if false, regrid only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution.
        :param v: If all voxels are the same, provides the initial voxel dimensions in h^-1 Mpc in dimensions (x, y, t), and auto-generates the voxel configuration array.
        :param osm_output: The relative path to save the osm file to.
        :param cosmology: The specific cosmology parameters in the form of a custom Cosmo object.
        :param save_dynamic_settings: If non-empty, save the dynamic settings to an .ini file given by the path entered.
        :param ref_time: An astropy.time.Time object stating the desired mid-observation time.
        :param ref_location: An astropy.coordinates.EarthLocation object stating the location of the telescope on Earth.
        :param observation_length: An astropy.time.TimeDelta object that gives the length of the observation.

        :return: The dynamically defined settings dictionary.
        """
        print("Initialising ...")

        # Configure d variable
        d = np.shape(values)

        # Set default voxel array according to v
        if voxels is None:
            print("Creating mock voxels ...")
            voxels = np.full((*d, 3), v, dtype=np.float64)

        # Transform datacube
        values, voxels, sigma_f, f_ref, regrid_flag = Regrid.transform_datacube_units(values=values, voxels=voxels, z_ref=z_ref, require_regrid=require_regrid, max_freq_res=max_freq_res, cosmology=cosmology)

        # STEP 7 - Regrid frequency-dimension data if needed
        if regrid_flag:
            values, voxels, sigma_f = Regrid.regrid_datacube(values=values, voxels=voxels, sigma_f=sigma_f, max_freq_res=max_freq_res)
        else:
            print("No regrid required!")

        # STEP 8 - Write data to OSM file
        Regrid.save_datacube_to_osm(values=values, voxels=voxels, sigma_f=sigma_f, f_ref=f_ref, phase_ref_point=phase_ref_point, osm_output=osm_output)

        # Output dynamic settings file
        dynamic_settings = Regrid.generate_dynamic_settings(
            values=values,
            voxels=voxels,
            f_ref=f_ref,
            phase_ref_point=phase_ref_point,
            ref_time=ref_time,
            ref_location=ref_location,
            observation_length=observation_length,
            save_dynamic_settings=save_dynamic_settings
        )

        return dynamic_settings
        
    @staticmethod
    def generate_osm_from_H5(file, phase_ref_point = RegridHelper.ZENITH_530, require_regrid = True, max_freq_res = 100e6, osm_output="regrid/osm_output/osm_output.osm", coeval=True, ref_time = RegridHelper.REF_TIME, ref_location = RegridHelper.SKA_REF_LOC, observation_length = RegridHelper.OBS_LEN_4HR, save_dynamic_settings = ""):
        """
        Combines both the convert_H5_to_csv and generate_osm_from_simulation functions.

        :param file: Location of the H5 file.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always regrid frequency bins, if false, regrid only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution in Hz.
        :param osm_output: The directory to output the osm file.
        :param ref_time: An astropy.time.Time object stating the desired mid-observation time.
        :param ref_location: An astropy.coordinates.EarthLocation object stating the location of the telescope on Earth.
        :param observation_length: An astropy.time.TimeDelta object that gives the length of the observation.
        :param save_dynamic_settings: If non-empty, save the dynamic settings to an .ini file given by the path entered.

        :return: The dynamically defined settings dictionary.
        """

        values, z_ref, vox, cosmology = Regrid.convert_H5_to_csv(file, coeval=coeval)

        if osm_output == "": osm_output = file.split('/')[-1][:-3] + "_osm.osm"

        return Regrid.generate_osm_from_simulation(values,
                z_ref=z_ref,
                require_regrid=require_regrid,
                max_freq_res=max_freq_res,
                v=vox, osm_output=osm_output,
                cosmology=cosmology,
                phase_ref_point=phase_ref_point,
                ref_time=ref_time,
                ref_location=ref_location,
                observation_length=observation_length,
                save_dynamic_settings=save_dynamic_settings
                )
    
    @staticmethod
    def convert_osm_file_to_arrays(osm_file, generate_dynamic_settings = True, phase_ref_point_override = None, ref_time = RegridHelper.REF_TIME, ref_location = RegridHelper.SKA_REF_LOC, observation_length = RegridHelper.OBS_LEN_4HR, save_dynamic_settings = "", d = None):
        """
        Reverse-engineer an osm file to retreive its values, voxels, sigma_f, f_ref, phase_ref_point, and dynamic settings.

        NB: The OSM file must be sorted according to the same order that it would've been constructed in i.e. frequency (descending), declination (ascending), right ascension (ascending).

        :param osm_file: The OSM file to analyse.
        :param generate_dynamic_settings: Whether or not to reverse-engineer the dynamic settings as well.
        :param phase_ref_point_override: If the phase refrence point is already known, override what the calculated phase refrence point would be.
        :param ref_time: An astropy.time.Time object stating the desired mid-observation time.
        :param ref_location: An astropy.coordinates.EarthLocation object stating the location of the telescope on Earth.
        :param observation_length: An astropy.time.TimeDelta object that gives the length of the observation.
        :param save_dynamic_settings: If non-empty, save the dynamic settings to an .ini file given by the path entered.
        :param d: The dimensions of the array, if none, the program will assume a cubic values array with side-lengths equal to the floor of the cube root of the number of entries in the OSM file.

        :return: The values, voxels, cumulative voxels, sigma_f, f_ref, and phase_ref_point contained in a dictionary. If generate_dynamic_settings is set to True, additionally return an updated dictionary of settings to provide to OSKAR.
        """

        df = pd.read_csv(osm_file, delimiter=" ", skiprows=3, index_col=False, names=["RA", "Dec", "Stokes I", "Q", "U", "V", "Freq0"])

        #values, voxels = None, cumulative_voxels = None, phase_ref_point = RegridHelper.ZENITH_530, f_ref = 200 * u.MHz
        output_data = {
            "values": None,
            "voxels": None,
            "cumulative_voxels": None,
            "phase_ref_point": None,
            "f_ref": None
        }

        # Get dimensions of box
        if d is None:
            d = tuple((np.ones(3, dtype=np.int32) * int(np.floor(np.cbrt(df.shape[0])))).tolist())

        # Extract values
        output_data["values"] = np.array(df["Stokes I"]).reshape(d)

        # Extract cumulative voxels
        output_data["cumulative_voxels"] = (np.array(df["RA"]).reshape(d), np.array(df["Dec"]).reshape(d), np.array(df["Freq0"]).reshape(d))

        # Create temporary values from the cumulative voxels
        temp_values = np.moveaxis(np.array(output_data["cumulative_voxels"]), 0, -1)

        # If the phase refrence point is given, override the calculation of the phase refrence point
        if phase_ref_point_override is None:
            output_data["phase_ref_point"] = SkyCoord(ra=Maths.denorm(Maths.circmean(temp_values[:,:,:,0]))*u.deg, dec=Maths.norm(Maths.circmean(temp_values[:,:,:,1], (-180, 180)))*u.deg, frame='icrs') 
        else:
            output_data["phase_ref_point"] = phase_ref_point_override

        # Calculate the frequency bin step
        step = np.mean(temp_values[:,:,0, 2] - temp_values[:,:,1, 2])/2

        # Calculate the refrence frequency
        output_data["f_ref"] = np.mean(temp_values[:,:,0, 2]) + step

        # Begin mutating temp_values

        # Shift values based on refrence points
        temp_values[:,:,:,0] = Maths.diff_sgn(temp_values[:,:,:,0], output_data["phase_ref_point"].ra.deg)
        temp_values[:,:,:,1] = Maths.diff_sgn(temp_values[:,:,:,1], output_data["phase_ref_point"].dec.deg)
        temp_values[:,:,:,2] = output_data["f_ref"] - temp_values[:,:,:,2]

        # Calculate the inverse cumulative sum
        temp_values[:,:,:,0] = np.diff(temp_values[:,:,:,0], prepend=0, axis=0)
        temp_values[:,:,:,1] = np.diff(temp_values[:,:,:,1], prepend=0, axis=1)
        temp_values[:,:,:,2] = np.diff(temp_values[:,:,:,2], prepend=-step, axis=2)

        # Extract voxels
        output_data["voxels"] = temp_values

        # Create temporary values from the cumulative voxels
        temp_values = np.moveaxis(np.array(output_data["cumulative_voxels"]), 0, -1)
        output_data["voxels"] = temp_values

        # Generate the dynamic settings
        if generate_dynamic_settings:
            dynamic_settings = Regrid.generate_dynamic_settings(
                values = output_data["values"],
                voxels = output_data["voxels"],
                cumulative_voxels = output_data["cumulative_voxels"],
                phase_ref_point = output_data["phase_ref_point"],
                f_ref = output_data["f_ref"] * u.Hz,
                ref_time=ref_time,
                ref_location=ref_location,
                observation_length=observation_length,
                save_dynamic_settings=save_dynamic_settings
                )

            return output_data, dynamic_settings
        
        else:
            return output_data


class BTAnalysisPipeline(object):
    """
    A broader class that combines all components of the individual components of the simulated IGM to simulated observation pipeline together.
    """

    @staticmethod
    def configure_oskar_settings(dynamic_settings = RegridHelper.DEFAULT_GENERAL_SETTINGS, interferometer_settings_override = "", imager_settings_override = "", save_ini=""):
        """
        Configure the settings files for the OSKAR interferometer and imager programs.

        :param osm_file: File path to the existing OSM file
        :param dynamic_settings: Dynamically generated settings from the regridding code.
        :param interferometer_settings_override: The file location of the OSKAR imager settings file. Leave blank if no override.
        :param imager_settings_override: The file location of the OSKAR interferometer settings template file. Leave blank if no override.
        :param use_imager: Whether or not to generate a dirty image with oskar_imager.
        :param save_ini: File path to save the compiled ini file. If blank, pass the settings only as a return.

        :return: The updated settings dictionary.
        """

        # Function to mutate the existing dynamic settings dictionary
        def mutate_settings(override_file, settings_dict):
            override_data = RegridHelper.read_settings_to_dictionary(override_file)

            return settings_dict | override_data | RegridHelper.PRIMARY_GENERAL_SETTINGS

        # Setup the interferometer ini file
        if interferometer_settings_override != "":
            dynamic_settings = mutate_settings(interferometer_settings_override, dynamic_settings)
            
        if imager_settings_override != "":
            dynamic_settings = mutate_settings(imager_settings_override, dynamic_settings)

        if save_ini != "":
            RegridHelper.save_settings_from_dictionary(save_ini, dynamic_settings)
        
        return dynamic_settings
    
    @staticmethod
    def split_general_settings(settings = ("", None), use_imager = True, save_file=True):
        """
        Split a settings dictionary into its components for the OSKAR interferometer and OSKAR imager.

        :param settings: A tuple containing a filepath and dictionary, the code will always prioritise using an evaluation the dictionary if available.
        :param use_imager: Whether to create a seperate imager file path and settings dictionary pair.
        :param save_file: If true, save to a file with a path and base name given by settings.

        :return: The two tuples containing the split settings (settings file location, settings dictionary)
        """

        # Declare default return values
        interf_settings_dict = None
        imager_settings_dict = None

        interf_settings_path = ""
        imager_settings_path = ""

        # If file is provided but not a dictionary, read the file
        if settings[1] is None and settings[0] != "":
            settings[1] = RegridHelper.read_settings_to_dictionary(settings[0])

        # Make sure that the settings exists before splitting
        if not settings[1] is None:
            interf_settings_dict = {k: settings[k] for k in RegridHelper.LEGAL_INTERFEROMETER_HEADINGS}

            if use_imager:
                imager_settings_dict = {k: settings[k] for k in RegridHelper.LEGAL_IMAGER_HEADINGS}

        # Make sure a file path has been provided to save the file to
        if settings[0] != "" and save_file:
            interf_settings_path = settings[0] + ".oskar_sim_interferometer.ini"

            RegridHelper.save_settings_from_dictionary(interf_settings_path, interf_settings_dict)

            if use_imager:
                imager_settings_path = settings[0] + ".oskar_imager.ini"
                
                RegridHelper.save_settings_from_dictionary(imager_settings_path, imager_settings_dict)

        return (interf_settings_path, interf_settings_dict), (imager_settings_path, imager_settings_dict)

    @staticmethod
    # FIXME: Read interferometer and imager settings files
    def run_oskar_on_osms(osm_file, interferometer_settings = ("", RegridHelper.DEFAULT_INTERFEROMETER_SETTINGS), imager_settings = ("", RegridHelper.DEFAULT_IMAGER_SETTINGS), fits_output="./fits_output.fits", oskar_exec=None, oskar_mode="python", use_imager=True):
        """
        Run oskar on each of the OSM sky models found in a fits directory, should already be formatted according to the output of the Regrid object.

        :param osm_file: Directory containing the OSM files to be imaged.
        :param interferometer_settings: A tuple containing the file path to the OSKAR interferometer settings file and the interferometer settings dictionary.
        :param imager_settings: A tuple containing the file path to the OSKAR imager settings file and the imager settings dictionary.
        :param fits_output: The directory and file to output the resultant FITS files.
        :param oskar_exec: The SIF file or binary file path containing the OSKAR programs.
        :param oskar_mode: How shall OSKAR be run? Options include: python, binary, command, singularity.
        :param use_imager: Whether or not to generate a dirty image with oskar_imager.
        """

        # TODO: Provide options for all four execution modes.
        # Create the output files
        subprocess.run(["mkdir","-p","BTA/oskar_output"], check=True)
        subprocess.run(["mkdir","-p","BTA/oskar_output/sim.ms"], check=True)

        print("Setting up OSKAR for "+osm_file)

        cwd = os.getcwd()+'/BTA'

        # Run OSKAR's interferometer simulation
        try:
            print("Running interferometer on "+osm_file)
            if oskar_mode == "singularity":
                subprocess.run(["singularity","exec","--nv","--bind",cwd,"--cleanenv","--home",cwd,oskar_exec,"oskar_sim_interferometer","sim_intif.ini"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
            elif oskar_mode == "binary":
                subprocess.run([oskar_exec+"/oskar_sim_interferometer","sim_intif.ini"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
        except subprocess.CalledProcessError as e:
            print(f"Command '{e.cmd}' returned non-zero exit status {e.returncode}.")
            print(f"Error output: {e.stderr.decode()}")

        # Run OSKAR's imager simulation
        if use_imager:
            try:
                print("Running imager on "+osm_file)
                if oskar_mode == "singularity":
                    subprocess.run(["singularity","exec","--nv","--bind",cwd,"--cleanenv","--home",cwd,oskar_exec,"oskar_imager","test_img.ini"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
                elif oskar_mode == "binary":
                    subprocess.run([oskar_exec+"/oskar_imager","test_img.ini"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
            except subprocess.CalledProcessError as e:
                print(f"Command '{e.cmd}' returned non-zero exit status {e.returncode}.")
                print(f"Error output: {e.stderr.decode()}")

        subprocess.run(["find",".","-name","'*.log'","-type","f","-delete"], check=True, cwd=cwd)
        subprocess.run(["cp","BTA/oskar_output/sim_image_I.fits",fits_output], check=True)

    @staticmethod
    def setup_bta_dir(h5_file, interferometer_settings_override="./regrid/test_intif_inis/test_img_gen.ini", imager_settings_override="./regrid/test_intif_inis/test_intif_gen.ini", oskar_telescope_model="./oskar_run_stage/telescope_model_AAstar", template=False):
        """
        Sets up the operating directory from which all anaysis will be done.

        :param h5_file: The location of the H5 file.
        :param cd_in: Whether to cd into the directory once finished or not.
        :param interferometer_settings_override: The file location of the OSKAR imager settings file. Leave blank if no override.
        :param imager_settings_override: The file location of the OSKAR interferometer settings template file. Leave blank if no override.
        :param oskar_telescope_model: The telescope model for OSKAR to use.
        :param template: If true, handle and return no h5 data.
        :return: The new location of the H5, imager ini, and interterometer template ini files, as well as the telescope model location. Also returns the PWD if cd_in is True.
        """

        cwd = os.getcwd()

        # 1. Create Directory
        subprocess.run(["mkdir","-p","BTA"], check=True)

        # 2. Move H5 file and INIs to directory
        if not template: subprocess.run(["cp",h5_file,"BTA/analysis.h5"], check=True)
        subprocess.run(["cp",interferometer_settings_override,"BTA/interferometer_override.ini"], check=True)
        subprocess.run(["cp",imager_settings_override,"BTA/imager_override.ini"], check=True)
        if not os.path.isdir("BTA/telescope_model"):
            subprocess.run(["cp","-r",oskar_telescope_model,"BTA/telescope_model"], check=True)

        return ("BTA/analysis.h5" if not template else ""), "BTA/interferometer_override.ini", "BTA/imager_override.ini", "BTA/telescope_model", cwd

    @staticmethod
    def clean_bta_dir(outdir, fits_output, clean=True):
        """
        Finishes and cleans up the mess created by the BTA class.

        :param outdir: The directory to save the finalised fits image, relative to parent dir.
        :param fits_cube: The name of the FITS output file.
        :param clean: Whether or not to remove the BTA directory. Only works if cd_out is true.
        """

        # Move datacube to parent directory
        subprocess.run(["mv", fits_output, outdir], check=True)

        if clean:
            subprocess.run(["rm","-rf","BTA"], check=True)

    @staticmethod
    def h5_box_to_datacube(file, phase_ref_point = RegridHelper.ZENITH_530, require_regrid = True, max_freq_res = 100e6, interferometer_settings_override = "./regrid/test_intif_inis/test_img_gen.ini", imager_settings_override = "./regrid/test_intif_inis/test_intif_gen.ini", outdir = ".", clean = True, oskar_exec = RegridHelper.OSKAR_SIF, oskar_mode="singularity", oskar_telescope_model = RegridHelper.TELESCOPE, template_preset = "", coeval = True, load_osm=False, ref_time = RegridHelper.REF_TIME, ref_location = RegridHelper.SKA_REF_LOC, observation_length = RegridHelper.OBS_LEN_4HR, use_imager = True):
        """
        Full pipeline function for transforming a H5 simulation box output into a FITS datacube.

        :param file: Location of the H5 file.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always regrid frequency bins, if false, regrid only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution in Hz.
        :param interferometer_settings_override: The file location of the OSKAR imager settings file.
        :param imager_settings_override: The file location of the OSKAR interferometer settings template file.
        :param outdir: The output location of the final FITS file.
        :param oskar_exec: The SIF file or location of compiled OSKAR binaries.
        :param oskar_mode: How shall OSKAR be run? Options include: python, binary, command, singularity.
        :param oskar_telescope_model: The telescope model for OSKAR to use.
        :param template_preset: Use a mock values array instead of a h5 file. Ignores any provided h5 file.
        :param coeval: If the H5 box is coeval or lightcone based.
        :param load_osm: If true load treat the file variable as if it were an OSM file.
        :param ref_time: An astropy.time.Time object stating the desired mid-observation time.
        :param ref_location: An astropy.coordinates.EarthLocation object stating the location of the telescope on Earth.
        :param observation_length: An astropy.time.TimeDelta object that gives the length of the observation.
        :param use_imager: Whether or not to generate a dirty image with oskar_imager.
        """

        # Expand paths
        if oskar_mode == "binary" or oskar_mode == "singularity":
            oskar_exec = RegridHelper.expand_path(oskar_exec)
        
        file = RegridHelper.expand_path(file)
        imager_settings_override = RegridHelper.expand_path(imager_settings_override)
        interferometer_settings_override = RegridHelper.expand_path(interferometer_settings_override)
        oskar_telescope_model = RegridHelper.expand_path(oskar_telescope_model)
        outdir = RegridHelper.expand_path(outdir)

        template_flag = (template_preset != "")
        
        # Set templates
        if load_osm: template_preset = file.split('/')[-1][:-4]

        print("Setting up BTA directory ...")
        h5_file, interf_override_ini, imager_override_ini, _, _ = BTAnalysisPipeline.setup_bta_dir(file,
            oskar_telescope_model=oskar_telescope_model,
            interferometer_settings_override=interferometer_settings_override,
            imager_settings_override=imager_settings_override,
            template=template_flag
            )

        # Create output file locations
        h5_id = file.split('/')[-1][:-3] if not template_flag else template_preset
        osm_output = RegridHelper.expand_path("BTA/" + h5_id + "_sky_model.osm")
        ini_output = RegridHelper.expand_path("BTA/" + h5_id + "_general_settings.ini")
        fits_output = RegridHelper.expand_path("BTA/" + h5_id + "_datacube.fits")

        # Set the default dynamic settings array
        dynamic_settings = RegridHelper.DEFAULT_GENERAL_SETTINGS
        
        # Run the OSM and Settings generators
        if not os.path.isfile(osm_output):
            if not load_osm:
                if template_flag:
                    # IF we want to generate a fresh osm file AND its from a template
                    print("Generating OSM files from template ...")
                    template_values = Regrid.mock_values(template_preset)
                    dynamic_settings = Regrid.generate_osm_from_simulation(
                        template_values,
                        phase_ref_point=phase_ref_point,
                        require_regrid=require_regrid,
                        max_freq_res=max_freq_res,
                        osm_output=osm_output,
                        ref_time=ref_time,
                        ref_location=ref_location,
                        observation_length=observation_length,
                        save_dynamic_settings=ini_output
                        )
                else:
                    # IF we want to generate a fresh osm file AND its from a provided h5 file
                    print("Generating OSM files from H5 ...")
                    dynamic_settings = Regrid.generate_osm_from_H5(
                        h5_file,
                        phase_ref_point=phase_ref_point,
                        require_regrid=require_regrid,
                        max_freq_res=max_freq_res,
                        osm_output=osm_output,
                        coeval=coeval,
                        ref_time=ref_time,
                        ref_location=ref_location,
                        observation_length=observation_length,
                        save_dynamic_settings=ini_output
                        )
            else:
                if template_flag:
                    # IF we want to skip generating the osm file AND a template osm has been specified
                    subprocess.run(["cp", RegridHelper.expand_path("~/.oskar/osm_templates/"+template_preset+"_sky_model.osm"), osm_output], check=True)
                    subprocess.run(["cp", RegridHelper.expand_path("~/.oskar/ini_templates/"+template_preset+"_general_settings.ini"), ini_output], check=True)

                    dynamic_settings = RegridHelper.read_settings_to_dictionary(ini_output)
                else:
                    # IF we want to skip generating the osm file AND a use a specified already-complete osm file
                    subprocess.run(["cp", file, osm_output], check=True)

                    _, dynamic_settings = Regrid.convert_osm_file_to_arrays(
                        osm_output,
                        generate_dynamic_settings=True,
                        phase_ref_point_override=phase_ref_point,
                        ref_time=ref_time,
                        ref_location=ref_location,
                        observation_length=observation_length,
                        save_dynamic_settings=ini_output
                        )

        # Configure the OSKAR settings
        dynamic_settings = BTAnalysisPipeline.configure_oskar_settings(
            dynamic_settings=dynamic_settings,
            interferometer_settings_override=interf_override_ini,
            imager_settings_override=imager_override_ini,
            save_ini=ini_output
            )

        # Split the current settings files and retreive specific files
        interferometer_settings, imager_settings = BTAnalysisPipeline.split_general_settings(settings=(ini_output, dynamic_settings), use_imager=use_imager, save_file=True)

        # Run OSKAR
        print("Running OSKAR ...")
        BTAnalysisPipeline.run_oskar_on_osms(
            osm_file=osm_output,
            interferometer_settings=interferometer_settings,
            imager_settings=imager_settings,
            fits_output=fits_output,
            oskar_exec=oskar_exec,
            oskar_mode=oskar_mode,
            use_imager=use_imager
            )

        # If clean is true remove all data relating to execution
        print("Cleaning up ...")
        BTAnalysisPipeline.clean_bta_dir(outdir=outdir, fits_output=fits_output, clean=clean)


# Testing stage
# pylint: disable=line-too-long

def load_defaults():
    """
    (Re)load default values
    """
    for template_presett in Regrid.TEMPLATE_PRESETS:
        if "coeval" in template_presett:
            template_value = Regrid.mock_values(template_presett, d=(400, 400, 400))
        else:
            template_value = Regrid.mock_values(template_presett, scale=20)

        Regrid.generate_osm_from_simulation(template_value, osm_output=RegridHelper.expand_path("~/.oskar/osm_templates/"+template_presett+"_sky_model.osm"), save_dynamic_settings=RegridHelper.expand_path("~/.oskar/ini_templates/"+template_presett+"_general_settings.ini"))

#load_defaults()

# FIXME: Test refactored code

#BTAnalysisPipeline.h5_box_to_datacube(None, template_preset="gaussian")

#BTAnalysisPipeline.h5_box_to_datacube("./regrid/osm_output/yuxiang1_zenith_osm", oskar_exec=RegridHelper.OSKAR_BIN, load_osm=True, oskar_mode="binary", oskar_telescope_model=RegridHelper.TELESCOPE)
