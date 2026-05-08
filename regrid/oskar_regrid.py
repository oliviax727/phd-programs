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

# TODO: Turn into pip project (later)

# TODO: Clean pylint errors
# pylint: disable=missing-function-docstring
# pylint: disable=unnecessary-lambda-assignment
# pylint: disable=invalid-name

# pylint: disable-next=unused-import
from matplotlib import pyplot as plt

# Regrid helper functions
class RegridHelper():
    """
        Helper functions and constants for the Regridding process.
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
        "general": {
            "app": "oskar_sim_interferometer"
        },
        "simulator": {
            "use_gpus": False
        },
        "observation" : {
            "num_time_steps": 24
        },
        "telescope": {
            "input_directory": "telescope_model",
            "apeture_array": {
                "element_pattern": {
                    "enable_numerical": "false"
                },
                "array_pattern": {
                    "element": {
                        "x_gain": 1.0,
                        "y_gain": 1.0,
                        "x_gain_error_time": 0.0015057,
                        "y_gain_error_time": 0.0015057,
                        "x_phase_error_fixed_deg": 0.0,
                        "y_phase_error_fixed_deg": 0.0,
                        "x_phase_error_time_deg": 0.0015057,
                        "y_phase_error_time_deg": 0.0015057
                    }
                }
            }
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
        "general": {
            "app": "oskar_imager"
        },
        "image": {
            "use_gpus": False,
            "channel_snapshots": "false",
            "input_vis_data": "output/sim.ms",
            "root_path": "output/sim_image"
        }
    }

    DEFAULT_GENERAL_SETTINGS = DEFAULT_IMAGER_SETTINGS | DEFAULT_INTERFEROMETER_SETTINGS | { "general": {} }

    # Angular distance calculation
    norm = lambda t: (t % 360 + 360) % 360 - 180
    delta = lambda a, b: RegridHelper.norm(a)-RegridHelper.norm(b)
    diff = lambda a, b: min(abs(RegridHelper.delta(a, b)),360-abs(RegridHelper.delta(a, b)))
    diff_sgn = lambda a, b: RegridHelper.delta(a, b) - 360 if RegridHelper.delta(a, b) > 180 else (RegridHelper.delta(a, b) + 360 if RegridHelper.delta(a, b) < -180 else RegridHelper.delta(a, b))

    # l, m, n to RA, Dec
    @staticmethod
    def lm_to_radec(l, m, phase_centre=ZERO_RADEC):
        d0 = phase_centre.dec.to_value(u.rad)
        a0 = phase_centre.ra.to_value(u.rad)
        n = np.sqrt(1-l**2-m**2)

        Dec = np.arcsin(m*np.cos(d0)+n*np.sin(d0))
        Ra = a0 + np.arctan(l/(n*np.cos(d0)-m*np.sin(d0)))

        return [Ra, Dec]

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

    @staticmethod
    def mock_values(preset, scale = 10, d = (100, 100, 100)):
        """
        Create an array of mock simulation values.

        :param preset: Mock brightness temperature array format. Options are {"flat", "random", "gaussian", "sinusoid"}
        :param scale: Define the Kelvin scale of the array (e.g. default value will create a uniform array of 10 K or a random array of 0 - 10 K).
        :param d: The size of the values datacube.
        :return: Mock brightness temperature values.
        """
        print("Creating mock brightness temperatures ...")

        values = np.zeros(d).astype(np.float64)

        for x in range(d[0]):
            for y in range(d[1]):
                for t in range(d[2]):
                    values[x, y, t] = np.sqrt((x-d[0]/2)**2 + (y-d[1]/2)**2)

        normal = lambda x, mean=0, var=1: np.exp(-(x-mean)**2/(2*var))/np.sqrt(2*np.pi*var)
        sinusoid = lambda x, f=1, ph=0: np.sin(2*np.pi*f*x+ph)

        # Set values array
        if preset == "random":
            # Completely random values
            values = np.random.rand(*d) * scale
        elif preset == "flat":
            # Completely uniform values
            values = np.ones(d) * scale
        elif preset == "gaussian":
            # Simple gaussian over sky in y-direction
            values = normal(values, var=d[2])
        elif preset == "sinusoid":
            # Simple sine wave over sky in x-direction
            values = sinusoid(values, f=5/d[2])
        elif preset == "point":
            # Point source with diffusion of 2-Sigma
            values = normal(values, var=2)

        # FIXME: Add OzStar H5 files as template files

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
    def transform_datacube_units(values, voxels, z_ref = 7, require_regrid = True, max_freq_res = 100 * u.MHz, v = (1, 1, 1), cosmology=Cosmo()):
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
        d = values.shape()

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

                    if (x == 0 and y == 0):
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
        d = values.shape()

        for x in range(d[0]):
            for y in range(d[1]):

                print("\rSpaxel # (", x, ",", y, ")", end="")

                # Set the values as being in the middle of each bin
                freq_values = np.cumsum(voxels[x, y, :, 2]) - voxels[x, y, :, 2]/2

                # Create interpolation B-spline
                bspline = misp(freq_values, values[x, y, :])
                dspline = misp(freq_values, sigma_f[x, y, :])

                # Check maximum frequency resolution
                if np.abs(freq_values[0], freq_values[-1])/d[2] > max_freq_res_hz:
                    # Resize d[2] to match maximum resolution
                    d[2] = np.abs(freq_values[0], freq_values[-1])/max_freq_res_hz

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

        # RA, Dec centering function
        # Add half the total and half the spaxel widths to centre the main point
        centering = (lambda x: (x - np.max(x, axis=(0, 1))/2 + np.min(x, axis=(0, 1))/2))

        rasum = centering(np.cumsum(voxels[:,:,:,0], axis=0))
        decsum = centering(np.cumsum(voxels[:,:,:,1], axis=1))
        freqsum = f_ref.to_value(u.Hz) - (np.cumsum(voxels[:,:,:,2], axis=2))

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
        d = values.shape()

        # Cumulative sums are more important than voxel bins now
        (RAs, Dcs, freqsum) = (None, None, None) # Keep Pylint Happy
        if cumulative_voxels is None and voxels is None:
            raise ValueError("Error: Either an array of voxels or cumulative voxes must be provided!")
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
                        if Dcs[x, y, t] >= 0: Decln = "+" + str(Decln)
                        else:                 Decln = "-" + str(Decln)

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
            warnings.warn("Warning: The provided observation time is too long! The object will not be in the sky for the entire duration of time.")

        return ref_time, obs_length_flag
    # pylint: enable=unused-argument

    @staticmethod
    def generate_dynamic_settings(values, voxels = None, cumulative_voxels = None, phase_ref_point = RegridHelper.ZENITH_530, f_ref = 200 * u.MHz, ref_time = RegridHelper.REF_TIME, ref_location = RegridHelper.SKA_REF_LOC, observation_length = RegridHelper.OBS_LEN_4HR):
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

        :return: The dynamically defined settings dictionary.
        """
        # Calculated settings: [observation] start_frequency_hz, num_channels, frequency_inc_hz, phase_centre_ra_deg, phase_centre_dec_deg, length, start_time_utc
        # Calculated settings: [image] fov_deg, size

        # SETUP
        # Cumulative sums are more important than voxel bins now
        (RAs, Dcs, freqsum) = (None, None, None) # Keep Pylint Happy
        if cumulative_voxels is None and voxels is None:
            raise ValueError("Error: Either an array of voxels or cumulative voxes must be provided!")
        elif voxels is None:
            (RAs, Dcs, freqsum) = cumulative_voxels
        elif cumulative_voxels is None:
            (RAs, Dcs, freqsum) = Regrid.calculate_cumulative_voxels(voxels=voxels, f_ref=f_ref, phase_ref_point=phase_ref_point)

        # Configure d variable
        d = values.shape()

        # Create deep copy of union/logical or settings set
        dynamic_settings = dict(RegridHelper.DEFAULT_INTERFEROMETER_SETTINGS)

        # BASIC CONFIG
        # Set starting frequency NB: The last channel has the lowest frequency!
        dynamic_settings['observation']['start_frequency_hz'] = np.mean(freqsum[:,:,-1])

        # Set number of channels and image size
        dynamic_settings['observation']['num_channels'] = d[2]
        dynamic_settings['image']['size'] = max(d[0], d[1])

        # Set the frequency increment
        dynamic_settings['observation']['start_frequency_hz'] = np.mean(voxels[:,:,:,2])

        # Set phase centre RA and Dec
        dynamic_settings['observation']['phase_centre_ra_deg'] = phase_ref_point.ra.deg
        dynamic_settings['observation']['phase_centre_dec_deg'] = phase_ref_point.dec.deg

        # COORDINATE CONFIG
        # Set image field of view and size
        ref_time, _ = Regrid.calculate_observation_time_from_date(phase_ref_point=phase_ref_point, ref_time=ref_time, ref_location=ref_location, observation_length=observation_length)

        # Calculate RA dimension
        rac = RegridHelper.diff(np.mean(RAs[-1,:,:]), np.mean(RAs[0,:,:]))

        # Calculate Dec dimension
        decc = RegridHelper.diff(np.mean(Dcs[-1,:,:]), np.mean(Dcs[0,:,:]))

        # Set the field of view
        dynamic_settings['image']['fov_deg'] = max(rac, decc)

        # Set the observation time and length
        dynamic_settings['observation']['start_time_utc'] = (ref_time - observation_length / 2).utc.value
        dynamic_settings['observation']['length'] = str(observation_length.to_value(format='datetime'))
        
        return dynamic_settings

    @staticmethod
    def generate_osm_from_simulation(values, voxels = None, z_ref = 7, phase_ref_point = RegridHelper.ZENITH_530, require_regrid = True, max_freq_res = 100 * u.MHz, v = (1, 1, 1), osm_output="regrid/osm_output/osm_output.osm", cosmology=Cosmo(), save_dynamic_settings = "", ref_time = RegridHelper.REF_TIME, ref_location = RegridHelper.SKA_REF_LOC, observation_length = RegridHelper.OBS_LEN_4HR):
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
        d = values.shape()

        # Set default voxel array according to v
        if voxels is None:
            print("Creating mock voxels ...")
            voxels = np.full((*d, 3), v, dtype=np.float64)

        # Transform datacube
        values, voxels, sigma_f, f_ref, regrid_flag = Regrid.transform_datacube_units(values=values, voxels=voxels, z_ref=z_ref, require_regrid=require_regrid, max_freq_res=max_freq_res, v=v, cosmology=cosmology)

        # STEP 7 - Regrid frequency-dimension data if needed
        if regrid_flag:
            values, voxels, sigma_f = Regrid.regrid_datacube(values=values, voxels=voxels, sigma_f=sigma_f, max_freq_res=max_freq_res)
        else:
            print("No regrid required!")

        # STEP 8 - Write data to OSM file
        Regrid.save_datacube_to_osm(values=values, voxels=voxels, sigma_f=sigma_f, f_ref=f_ref, phase_ref_point=phase_ref_point, osm_output=osm_output)

        # Output dynamic settings file
        dynamic_settings = Regrid.generate_dynamic_settings(values=values, voxels=voxels, f_ref=f_ref, phase_ref_point=phase_ref_point, ref_time=ref_time, ref_location=ref_location, observation_length=observation_length)

        if save_dynamic_settings != "":
            settings_path = RegridHelper.expand_path(save_dynamic_settings)

            config = cfp.ConfigParser()
            config.read_dict(dynamic_settings)

            with open(save_dynamic_settings, 'w', encoding='utf-8') as settings_file:
                config.write(settings_file)

            print("Saved dynamic and default settings to ini file: "+settings_path)

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
    
    # FIXME: Reverse-read OSM
    @staticmethod
    def convert_osm_file_to_arrays(osm_file, generate_dynamic_settings = True):
        """
        Reverse-engineer an osm file to retreive its values, voxels, sigma_f, f_ref, phase_ref_point, and dynamic settings.

        :param osm_file: The OSM file to analyse.
        :param generate_dynamic_settings: Whether or not to reverse-engineer the dynamic settings as well

        :return: The values, voxels, sigma_f, f_ref, and phase_ref_point contained in a tuple. If generate_dynamic_settings is set to True, additionally return an updated dictionary of settings to provide to OSKAR.
        """

        #df = pd.read_csv(osm_file, delimiter=" ", skiprows=3, index_col=False, names=["RA", "Dec", "Stokes I", "Q", "U", "V", "Freq0"])

        print(osm_file)

        output_data = (None, None, None, None, None)

        if generate_dynamic_settings:
            dynamic_settings = RegridHelper.DEFAULT_GENERAL_SETTINGS

            return output_data, dynamic_settings
        else:
            return output_data


class BTAnalysisPipeline(object):
    """
    A broader class that combines all components of the individual components of the simulated IGM to simulated observation pipeline together.
    """
    
    # FIXME: Implement config parser
    @staticmethod
    def configure_oskar_settings(osm_file, dynamic_settings = RegridHelper.DEFAULT_GENERAL_SETTINGS, interferometer_settings_override = "", imager_settings_override = "", use_imager=True, save_ini=""):
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

        #_, osm_settings = Regrid.convert_osm_file_to_arrays(osm_file)

        # Setup the interferometer ini file
        # FIXME: Write to new interferometer file and copy over settings
        #with open('sim_intif.ini', 'w', encoding='utf-8'):
        
            # Set sky model location
            #ofname = "oskar_sky_model/file="+osm_file

            # Set frequency bin
            #freq = osm.split("_")[-1][:-7]+"e6"
            #ofname = r"start_frequency_hz="+freq
            #BTAnalysisPipeline.find_replace_line("BTA/test_intif.ini", "freqset", ofname)

        if use_imager:
            print("imager")
            # Setup the imager ini file
            # FIXME: Write to new imager file and copy over settings
            #with open('sim_img.ini', 'w', encoding='utf-8'):
                # Set the pixel resolution
                #ofname = "size="+str(size) 
                #BTAnalysisPipeline.find_replace_line("BTA/test_img.ini", "sizeset", ofname)

                # Set the FOV
                #ofname = "fov_deg="+str(fov)
                #BTAnalysisPipeline.find_replace_line("BTA/test_img.ini", "fovset", ofname)

        if save_ini:
            print(0)
        
        return {}
    
    @staticmethod
    def split_general_settings(dynamic_settings = None, settings_file = "", use_imager = True):
        """
        Split a settings dictionary into its components for the OSKAR interferometer and OSKAR imager.

        :param dynamic_settings: Dynamically generated settings from the regridding code. If none do not mutate.
        :param settings_file: A file containing the ini settings. If empty don't split the settings files.
        :param use_imager: Whether to create a seperate imager file path and settings dictionary pair.

        :return: The two tuples containing the split settings (settings file location, settings dictionary)
        """

        interf_settings_dict = None
        imager_settings_dict = None

        interf_settings_path = ""
        imager_settings_path = ""

        if not dynamic_settings is None:
            print(0)

        if settings_file != "":
            print(0)

        if use_imager:
            return (interf_settings_path, interf_settings_dict), (imager_settings_path, imager_settings_dict)
        else:
            return (interf_settings_path, interf_settings_dict), (None, "")


    @staticmethod
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
        :param require_regrid: If true then always regrid frequency bins, if false,
        regrid only when max frequency resolution is met.
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
                    dynamic_settings = Regrid.generate_osm_from_simulation(template_values, osm_output=osm_output, save_dynamic_settings=ini_output)
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

                    config = cfp.ConfigParser()
                    config.read_file(open(ini_output, encoding='utf-8'))

                    dynamic_settings = {s:dict(config.items(s)) for s in config.sections()}
                else:
                    # IF we want to skip generating the osm file AND a use a specified already-complete osm file
                    subprocess.run(["cp", file, osm_output], check=True)

        # Configure the OSKAR settings
        dynamic_settings = BTAnalysisPipeline.configure_oskar_settings(
            osm_file=osm_output,
            dynamic_settings=dynamic_settings,
            interferometer_settings_override=interf_override_ini,
            imager_settings_override=imager_override_ini,
            use_imager=use_imager,
            save_ini=ini_output
            )

        # Split the current settings files and retreive specific files
        interferometer_settings, imager_settings = BTAnalysisPipeline.split_general_settings(dynamic_settings=dynamic_settings, settings_file=ini_output, use_imager=use_imager)

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

#Regrid.generate_osm_from_H5("./regrid/yuxiang_bts/yuxiang1.h5", osm_output="./regrid/osm_output/yuxiang1_non_uniform_zenith_osm", coeval=True, require_regrid=False)
#Regrid.generate_osm_from_H5("./regrid/yuxiang_bts/yuxiang1.h5", osm_output="./regrid/osm_output/yuxiang1_00_osm", coeval=True, phase_ref_point=RegridHelper.ZERO_RADEC)
#Regrid.generate_osm_from_H5("./regrid/yuxiang_bts/yuxiang1.h5", osm_output="./regrid/osm_output/yuxiang1_zenith_osm", coeval=True)

#for template_preset in ["gaussian", "point", "random", "flat", "sinusoid", "point"]:
#    template_value = Regrid.mock_values(template_preset, scale=20)
#    Regrid.generate_osm_from_simulation(template_value, osm_output="./regrid/osm_output/"+template_preset+"_zenith_osm")

#Collator.collate_fits("./regrid/test_output/yuxiang1_fits", "./regrid/test_output")
#Collator.collate_fits("./regrid/test_output/yuxiangbad_fits", "./regrid/test_output")

#BTAnalysisPipeline.h5_box_to_datacube(None, template_preset="gaussian")

BTAnalysisPipeline.h5_box_to_datacube("./regrid/osm_output/yuxiang1_zenith_osm", oskar_exec=RegridHelper.OSKAR_BIN, load_osm=True, oskar_mode="binary", oskar_telescope_model=RegridHelper.TELESCOPE)
