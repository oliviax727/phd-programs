"""
The oskar_helpers module contains a series of helper functions that handle the basic mathematical and computer methods involved in the oskareor.data module.
"""

# System imports
import os
import configparser as cfp
from typing import Iterable

# Mathematics and calculations
import astropy.constants as c
import astropy.units as u
import numpy as np

# Astropy extras
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.coordinates.errors import UnknownSiteException as USE
from astropy.cosmology import FlatLambdaCDM as fmodel
from astropy.cosmology import z_at_value as getz
from astropy.time import Time, TimeDelta
from astropy.units import Quantity

# Data handling and statistics
from scipy.stats import circmean as circmean_radians

class SKAMath():
    """
    Mathematical helper functions for various calculations used throughout the document. All angle units are in radians unless if otherwise specified.
    """

    # Retreive SKA Location from astropy server, if it fails, return the default Earth Location
    @staticmethod
    def get_ska_location() -> EarthLocation:
        """
        Gets the EarthLocation object for the SKA (or MWA).
        """

        ska_location = None

        try:
            ska_location = EarthLocation.of_site("SKA-LOW")
        except USE:
            try:
                ska_location = EarthLocation.of_site("MWA")
            except USE:
                ska_location = EarthLocation.from_geodetic(lon=-117*u.deg, lat=-27*u.deg)
        
        return ska_location

    # Astropy time/location/coordinate locations
    SKA_REF_LOC: EarthLocation  = get_ska_location()
    OBS_LEN_4HR: TimeDelta      = TimeDelta(4 * u.hr)
    REF_TIME: Time              = Time(val="2025-03-03T05:30:00.00", format='isot', scale='utc')
    ZENITH_530: SkyCoord        = SkyCoord(ra=0*u.deg, dec=-27*u.deg, frame='icrs') # SKA-Low Zenith at 5:30 am 2025-03-03
    ZERO_RADEC: SkyCoord        = SkyCoord(ra=0*u.deg, dec=0*u.deg, frame='icrs') # Centre RA/dec

    # Astropy unit constants
    FREQ_21CM: Quantity         = 1.420405751768e9 * u.Hz
    WLEN_21CM: Quantity         = (c.c / FREQ_21CM).to(u.cm)
    SIGMA_F: Quantity           = (np.sqrt(c.k_B / (1.008 * c.u * WLEN_21CM**2))).to(u.Hz*u.K**-0.5)

    # Unitless quantities
    SIGMA_F_FLOAT: float        = SIGMA_F.to_value(u.Hz * u.K ** (-1/2))

    # Angular distance calculations in degrees
    @staticmethod
    def normalise_angle(t: float) -> float:
        """ Calculate the normalised angle """
        return (t % 360 + 180) % 360 - 180
    
    normalise = normalise_angle
    norm = normalise_angle

    @staticmethod
    def denormalise_angle(t: float) -> float:
        """ Denormalise an angle i.e. convert to [0, 360). """
        return (t % 360 + 360) % 360
    
    denormalise = denormalise_angle
    denorm = denormalise_angle

    @staticmethod
    def angle_delta(a: float, b: float) -> float:
        """ Calculate the unnormalised difference between two angles. """
        return SKAMath.norm(a)-SKAMath.norm(b)

    delta = angle_delta

    @staticmethod
    def angle_difference(a: float, b: float) -> float:
        """ Calculate the true difference between two angles. """
        return min(abs(SKAMath.delta(a, b)),360-abs(SKAMath.delta(a, b)))

    angle_difference = np.vectorize(angle_difference)
    diff = angle_difference

    @staticmethod
    def signed_angle_difference(a: float, b: float) -> float:
        """ Calculate the signed difference between two angles. """
        return SKAMath.norm(SKAMath.delta(a, b) - 360) if (SKAMath.delta(a, b)) > 180 else ((SKAMath.delta(a, b) + 360) if SKAMath.delta(a, b) < -180 else SKAMath.delta(a, b))

    signed_angle_difference = np.vectorize(signed_angle_difference)
    diff_sgn = signed_angle_difference

    # Scipy circmean function for degrees
    @staticmethod
    def circmean_deg(angles: float, bounds: tuple[float] = (0, 360)) -> float:
        """ Calculate the mean using modular arithmetic, a modified version of scipy.stats.circmean. """
        return np.rad2deg(circmean_radians(np.deg2rad(angles), high=np.deg2rad(bounds[1]), low=np.deg2rad(bounds[0])))

    circmean = circmean_deg

    # Normal, sinusoid, and sinc functions for convenience
    @staticmethod
    def normal(x: float, mean: float = 0, var: float = 1, amp: float = np.sqrt(2*np.pi)) -> float:
        """ A simple normal distribution function. """
        return amp * np.exp(-(x-mean)**2/(2*var))/np.sqrt(2*np.pi*var)
    
    gaussian = normal
    
    @staticmethod
    def sinusoid(x: float, f: float = 1, ph: float = 0, amp: float = 1) -> float:
        """ A simple sinusoid function. """
        return amp * np.cos(2*np.pi*f*x+ph)

    @staticmethod
    def sinc(x: float, f: float = 1, ph: float = 0, amp: float = 1) -> float:
        """ A simple sinc function. """
        return amp * np.nan_to_num(np.sin(2*np.pi*f*x+ph)/(2*np.pi*f*x+ph), nan=1, posinf=1, neginf=1)
    
    # Convert FWHM to Variance and vice versa
    @staticmethod
    def full_width_at_half_maximum(x: float) -> float:
        """ Calculate the FWHM from the standard deviation. Distribution is assumed as gaussian. """
        return x * 2 * np.sqrt(2*np.log(2))
    
    FWHM = full_width_at_half_maximum

    @staticmethod
    def standard_deviation(x: float) -> float:
        """ Calculate the standard deviation from the FWHM. Distribution is assumed as gaussian. """
        return x / (2 * np.sqrt(2*np.log(2)))
    
    STDEV = standard_deviation

    # l, m, n to RA, dec
    @staticmethod
    def lm_to_radec(l: float, m: float, phase_centre: SkyCoord = ZERO_RADEC) -> list[float]:
        """ Convert the l, m plane coordinates to RA and dec. """
        d0 = phase_centre.dec.to_value(u.rad)
        a0 = phase_centre.ra.to_value(u.rad)
        n = np.sqrt(1-l**2-m**2)

        dec = np.arcsin(m*np.cos(d0)+n*np.sin(d0))
        ra = a0 + np.arctan(l/(n*np.cos(d0)-m*np.sin(d0)))

        return (ra, dec)

class EoRCosmology():
    """
    Defines a specific cosmological model for other classes to refer to. Assumes a flat universe.

    :param h0 (Quantity): The Hubble constant.
    :param omega_m_0 (float): The dimensionless matter density.
    :param omega_b_0 (float): The dimensionless baryonic matter density.
    :param cosmo (EoRCosmology): The cosmological ΛCDM model.
    """

    def __init__(self, h0: Quantity = 100 * u.km / u.s / u.Mpc, omega_m_0: float = 0.31, omega_b_0: float = 0.048):
        self.h0        = h0 # Set Hubble Constant to 100 h, with h being dimensionless hubble parameter
        self.omega_m_0 = omega_m_0
        self.omega_b_0 = omega_b_0
        self.cosmo     = fmodel(H0=h0, Om0=omega_m_0, Ob0=omega_b_0) # Flat ΛCDM model

    # Redshift to comoving distance
    def z_to_dz(self, z: float) -> Quantity:
        """ Convert redshift to comoving distance. """
        return self.cosmo.comoving_distance(z)

    # Comoving distance to redshift
    def dz_to_z(self, dz: Quantity) -> float:
        """ Convert comoving distance to redshift. """
        return getz(self.cosmo.comoving_distance, dz).value

    # Redshift to frequency in GHz
    @staticmethod
    def z_to_f(z: Quantity, f0: Quantity = 1.420e9 * u.Hz) -> Quantity:
        """ Convert redshift to frequency. By default it assumes that the start frequency is the 21 cm line. """
        return f0 / (z + 1)
    
    # Frequency in GHz to redshift
    @staticmethod
    def f_to_z(f: Quantity, f0: Quantity = 1.420e9 * u.Hz) -> Quantity:
        """ Convert frequency to redshift. By default it assumes that the start frequency is the 21 cm line. """
        return (f0 / f) - 1

class OSKARFileConfig:
    """
    The OSKARFileConfig class contains helper functions specifically relating to directory and file management.
    """

    @staticmethod
    def expand_path(path: str) -> str:
        """ Expands a filepath to produce an absolute path. """
        if path != "":
            return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
        
        raise ValueError("Cannot expand an empty path string.")

    @staticmethod
    def dir_list_sorted(dir_: str) -> str:
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
    def find_replace_line(file_name: str, find_line: str, replace_line: str):
        """
        Replace a given line in a settings.ini file given the line is equal to a special string.
        
        :param file (str): The file to perform the find-and-replace.
        :param find_line (str): The special keyword to trigger a replace.
        :param replace_line (str): The line to replace the preset.
        """

        with open(OSKARFileConfig.expand_path(file_name), 'r', encoding='utf-8') as file:
            lines = file.readlines()

        with open(OSKARFileConfig.expand_path(file_name), 'w', encoding='utf-8') as file:
            for line in lines:
                if line.startswith(find_line):
                    file.write(replace_line+"\n")
                else:
                    file.write(line)

    @staticmethod
    def read_settings_to_dictionary(file: str) -> dict:
        """
        Convert a settings ini file into a dictionary using configparser.

        :param file (str): The specific ini file to read.

        :return settings_dict (dict): The settings in the form of a dictionary.
        """

        config = cfp.ConfigParser()
        with open(OSKARFileConfig.expand_path(file), encoding='utf-8') as filestream:
            config.read_file(filestream)

        return { s: dict(config.items(s)) for s in config.sections() }
    
    @staticmethod
    def save_settings_from_dictionary(file: str, settings_dict: dict):
        """
        Save a dictionary contents as an ini file using configparser.

        :param file (str): The specific ini file to save the data to.
        :param settings_dict (dict): The settings in the form of a dictionary.
        """

        config = cfp.ConfigParser()
        config.read_dict(settings_dict)

        with open(OSKARFileConfig.expand_path(file), 'w', encoding='utf-8') as settings_file:
            config.write(settings_file, space_around_delimiters=False)

    @staticmethod
    def recase_iterable(iterable: Iterable[str], lower: bool = True) -> Iterable[str]:
        """ Makes every string in an iterable lower case then turns it back into the iterable object. """

        iter_type = type(iterable)
        map_obj = map(lambda s: s.lower() if lower else s.upper(), iterable)

        return iter_type(map_obj)

# TODO: Compile the ~/oskareor.data directory with the ./oskareor directory into one github project
