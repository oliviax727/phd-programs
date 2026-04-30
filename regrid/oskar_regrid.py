#!.venv/bin/python

# Import the stuffs
import os
import subprocess
import astropy.constants as c
import astropy.units as u
import h5py
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM as fmodel
from astropy.cosmology import z_at_value as getz
from scipy.interpolate import make_interp_spline as misp
import configparser as cfp

# FIXME: Turn into pip project (later)

# pylint: disable-next=unused-import
from matplotlib import pyplot as plt

# Regrid helper functions
class Regrid_Helper(object):
    """
        Helper functions and constants for the Regridding process.
    """

    # Define Constants
    ZENITH_530  = SkyCoord(ra=0*u.deg, dec=-27*u.deg, frame='icrs') # SKA-Low Zenith at 5:30 am 2025-03-03
    ZERO_RADEC  = SkyCoord(ra=0*u.deg, dec=0*u.deg, frame='icrs') # Centre RA/Dec
    OSKAR_SIF   = "~/.oskar/OSKAR-2.12.2-Python3.sif"
    OSKAR_BIN   = "~/.oskar/bin/"
    TELESCOPE   = "~/.oskar/SKA-Low_telescope_models/SKA-Low_AAstar_original_rigid-rotation.tm"
    SIGMA_F     = (np.sqrt(c.k_B / (1.008 * c.u * (21.106 * u.cm)**2))).to(u.Hz*u.K**-0.5).value

    # Define default settings
    # FIXME: Edit default settings
    DEFAULT_INI = {
        "simulator": {
            "use_gpus": False
        },
        "observation" : {
            "num_channels": 1,
            "start_frequency_hz": 177.408e6,
            "frequency_inc_hz": 20e6,
            "phase_centre_ra_deg": 0.0,
            "phase_centre_dec_deg": -27.0,
            "num_time_steps": 8,
            "start_time_utc": "2025-03-03 05:30:00.00",
            "length": "00:02:00.000"
        },
        "telescope": {
            "input_directory": "../oskar_run_stage/telescope_model_AAstar"
        },
        "interferometer": {
            "oskar_vis_filename": "oskar_output/vis.vis",
            "ms_filename": "oskar_output/sim.ms",
            "channel_bandwidth_hz": 1e6,
            "time_average_sec": 10.0
        },
        "sky": {}
    }

    # Angular distance calculation
    norm = lambda t: (t % 360 + 360) % 360 - 180
    delta = lambda a, b: Regrid_Helper.norm(a)-Regrid_Helper.norm(b)
    diff = lambda a, b: min(abs(Regrid_Helper.delta(a, b)),360-abs(Regrid_Helper.delta(a, b)))
    diff_sgn = lambda a, b: Regrid_Helper.delta(a, b) - 360 if Regrid_Helper.delta(a, b) > 180 else (Regrid_Helper.delta(a, b) + 360 if Regrid_Helper.delta(a, b) < -180 else Regrid_Helper.delta(a, b))

    # l, m, n to RA, Dec
    @staticmethod
    def lm_to_radec(l, m, phase_centre=Regrid_Helper.ZERO_RADEC):
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

class Cosmo(object):
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

class Regrid(object):
    """
    The regridding class contains functions relating to translating simulation data to OSKAR output data.
    """

    @staticmethod
    def brightness_temperature_to_flux(Tb, fxy, dθ, dφ):
        # Calculate pixelated luminosity
        Fv = 2 * c.k_B.value * fxy**2 * Tb * (dθ * dφ) / c.c.value ** 2
        Fv = Fv * 1e26 # Converts to Jansky
        return Fv
    
    @staticmethod
    def brightness_temperature_to_linewidth(Tb):
        # Calculate linewidth in Hz
        return (Tb ** 0.5) * Regrid_Helper.SIGMA_F

    @staticmethod
    def mock_values(preset, scale = 10, d = (100, 100, 100)):
        """
        Create an array of mock simulation values.

        :param preset: Mock brightness temperature array format. Options are {"flat", "random", "gaussian", "sinusoid"}
        :param scale: Define the Kelvin scale of the array (e.g. default value will create a uniform array of 10 K or a random array of 0 - 10 K).
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

    # FIXME: Seperate generate_osm_from_simulation into subfunctions
    @staticmethod
    def generate_osm_from_simulation(values, voxels = None, d = (100, 100, 100), z_ref = 7, phase_ref_point = Regrid_Helper.ZENITH_530, require_regrid = True, max_freq_res = 100e6, v = (1, 1, 1), osm_output="osm_output", cosmology=Cosmo()):
        """
        Generate a set of .osm files for an OSKAR sky model based on a Mpc**3 simulation output.

        :param values: The simulation datacube, must have shape d.
        :param voxels: An array describing a series of voxel dimensions corresponding to each simulation datacube voxel element.
        :param d: Number of voxels in simulation in dimensions (x, y, t).
        :param z_ref: Refrence redshift, the ending redshift of the simulation.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always regrid frequency bins, if false, regrid only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution in Hz.
        :param v: If all voxels are the same, provides the initial voxel dimensions in h^-1 Mpc in dimensions (x, y, t), and auto-generates the voxel configuration array.
        :param osm_output: The relative path to save the osm file to.
        :param cosmology: The specific cosmology parameters in the form of a custom Cosmo object.
        """
        print("Initialising ...")

        # Calculate refrence comoving dist.
        Dz_ref = cosmology.z_to_Dz(z_ref)
        f_ref = cosmology.z_to_f(z_ref)

        # Set voxels array
        if voxels == None:
            print("Creating mock voxels ...")
            voxels = np.full((*d, 3), v, dtype=np.float64)

        # Set regrid flag
        regrid_flag = require_regrid

        # Set cosmological redshift parameters
        Dz = Dz_ref
        z_prev = z_ref
        fq = f_ref.to_value(u.Hz)
        Dz_val = Dz_ref.to_value(u.Mpc)

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
                    dθ = np.arctan(dx/Dz_pix)
                    dφ = np.arctan(dy/Dz_pix)

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
                        if df > max_freq_res:
                            regrid_flag = regrid_flag or True

                    else:
                        df = voxels[0, 0, t, 2]

                    # STEP 4 - Convert brightness temperature to pixel flux
                    fxy = fq + df/2
                      
                    # Calculate pixelated luminosity
                    Fv = Regrid.brightness_temperature_to_flux(Tb=Tb, fxy=fxy, dθ=dθ, dφ=dφ)
                
                    # Save values
                    voxels[x, y, t, 0] = dθ
                    voxels[x, y, t, 1] = dφ
                    voxels[x, y, t, 2] = df
                    values[x, y, t] = Fv
            
            # Increment distance and frequency parameters
            Dz = Dz + dt*u.Mpc # Increment the value of Dz by voxel dimension
            z_prev = z_top # Top z in current box = bottom z in next box
            Dz_val = Dz_val + dt # Increment the value of Dz by voxel dimension
            fq = fq + df # Increment cumulative frequency

            print("\rTime step #", t, end="")

        print("\nTransforming complete.")

        # STEP 7 - Regrid frequency-dimension data if needed
        if regrid_flag:
            print("Performing regrid ...")

            for x in range(d[0]):
                for y in range(d[1]):

                    print("\rSpaxel # (", x, ",", y, ")", end="")

                    # Set the values as being in the middle of each bin
                    freq_values = np.cumsum(voxels[x, y, :, 2]) - voxels[x, y, :, 2]/2

                    # Create interpolation B-spline
                    bspline = misp(freq_values, values[x, y, :])
                    dspline = misp(freq_values, sigma_f[x, y, :])

                    # Perform regridding
                    new_freq, freq_bandw = np.linspace(freq_values[0], freq_values[-1], d[2], retstep=True)
                    new_flux = np.clip(bspline(new_freq), 0, None)
                    new_sigf = np.clip(dspline(new_freq), 0, None)

                    # Create array of uniform bin sizes
                    freq_bins = np.ones(d[2]) * freq_bandw

                    # Save variables
                    voxels[x, y, :, 2] = freq_bins
                    values[x, y, :] = new_flux
                    sigma_f[x, y, :] = new_sigf

            print("\nRegrid complete.")

        else:
            print("No regrid required!")
        
        # STEP 8 - Write data to OSM file
        print("Configuring datacube for OSKAR file format ...")

        # RA, Dec centering function
        # Add half the total and half the spaxel widths to centre the main point
        centering = (lambda x: (x - np.max(x, axis=(0, 1))/2 + np.min(x, axis=(0, 1))/2))

        # Cumulative sums are more important than voxel bins now
        rasum = centering(np.cumsum(voxels[:,:,:,0], axis=0))
        decsum = centering(np.cumsum(voxels[:,:,:,1], axis=1))
        freqsum = f_ref.to_value(u.Hz) - (np.cumsum(voxels[:,:,:,2], axis=2))

        # Calculate phase centre offsets
        source_pos = phase_ref_point.spherical_offsets_by(rasum * u.rad, decsum * u.rad)
        RAs = source_pos.ra.to_value(u.deg)
        Dcs = source_pos.dec.to_value(u.deg)

        # Record data to file
        print("Recording data to .osm file")

        with open(osm_output, 'w', encoding='utf8') as osm:
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
                        RAscn = np.char.zfill(np.format_float_positional(RAs[x, y, t], 6, False), 10)
                        Decln = np.char.zfill(np.format_float_positional(np.abs(Dcs[x, y, t]), 6, False), 9)
                        value = np.format_float_scientific(values[x, y, t], 4, False)
                        freq0 = np.format_float_positional(freqsum[x, y, t] / 1e6, 3, False)
                        linew = np.format_float_scientific(sigma_f[x, y, t], 4, False)

                        # Add +/- value to Declinations
                        if Dcs[x, y, y] >= 0: Decln = "+" + str(Decln)
                        else:                 Decln = "-" + str(Decln)

                        # Write to OSM
                        osm.write(
                            str(RAscn)    + " " + # Right Ascension
                            Decln         + " " + # Declination
                            str(value)    + " " + # Intensity
                            "0.0 0.0 0.0" + " " + # Redundant Stokes Parameters
                            str(freq0)  + "e6 " + # Point source frequency
                            str(linew)    + " " + # Spectral profile linewidth
                            "\n"
                        )

                    print("\rSpaxel # (", x, ",", y, ")", end="")

        print("\nProcess complete, data saved to "+osm_output)

    
    @staticmethod
    def generate_osm_from_H5(file, phase_ref_point = Regrid_Helper.ZENITH_530, require_regrid = True, max_freq_res = 100e6, osm_output="", coeval=True):
        """
        Combines both the convert_H5_to_csv and generate_osm_from_simulation functions.

        :param file: Location of the H5 file.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always regrid frequency bins, if false, regrid only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution in Hz.
        :param osm_output: The directory to output the osm file.
        """

        values, dim, z_ref, vox, cosmology = Regrid.convert_H5_to_csv(file, coeval=coeval)

        if osm_output == "": osm_output = file.split('/')[-1][:-3] + "_osm"

        Regrid.generate_osm_from_simulation(values, d=dim, z_ref=z_ref, require_regrid=require_regrid, max_freq_res=max_freq_res, v=vox, osm_output=osm_output, cosmology=cosmology, phase_ref_point=phase_ref_point)

class BTAnalysisPipeline(object):
    """
    A broader class that combines all components of the individual components of the simulated IGM to simulated observation pipeline together.
    """

    @staticmethod
    def get_osm_sky_dimensions(osm_file):
        """
        Get the pixel count for a FITS file with it's corresponding OSKAR sky model.

        :param osm_file: The OSM file to analyse.
        :return: The ideal image size and fov
        """

        df = pd.read_csv(osm_file, delimiter=" ", skiprows=3, index_col=False, names=["RA", "Dec", "Stokes I", "Q", "U", "V", "Freq0"])

        print(osm_file)

        # Calculate RA dimension
        RAc = Regrid_Helper.diff(np.array(df['RA'])[-1], np.array(df['RA'])[0])

        # Calculate Dec dimension
        Decc = Regrid_Helper.diff(np.array(df['Dec'])[-1], np.array(df['Dec'])[0])

        # Calculate FOV
        fov = max(RAc, Decc)

        # Get number of voxels on a side
        n = np.sqrt(len(np.array(df['RA'])))

        return n, fov
    
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
    def run_oskar_on_osms(osm_file, imager_template_ini = "./regrid/test_intif_inis/test_img_gen.ini", interferometer_template_ini = "./test_intif_inis/test_intif_gen.ini", fits_output="./fits_output.fits", oskar_exec=Regrid_Helper.OSKAR_SIF, oskar_mode="singularity"):
        """
        Run oskar on each of the OSM sky models found in a fits directory, should already be formatted according to the output of the Regrid object.

        :param osm_dir: Directory containing the OSM files to be imaged.
        :param imager_template_ini: The file location of the OSKAR imager settings file.
        :param interferometer_template_ini: The file location of the OSKAR interferometer settings template file. The ini file must contain a line under [observation] with text "fset" in liu of the start_frequency_hz setting, and a line under [sky] with text "preset" in liu of the oskar_sky_model/file setting.
        :param fits_output: The directory and file to output the resultant FITS files.
        :param oskar_exec: The SIF file containing the OSKAR program. OSKAR must be run from singularity.
        :param oskar_mode: How shall OSKAR be run? Options include: python, binary, command, singularity.
        """
        # FIXME: Provide options for all four execution modes.
        # Create the output files
        subprocess.run(["mkdir","-p","BTA/oskar_output"], check=True)
        subprocess.run(["mkdir","-p","BTA/oskar_output/sim.ms"], check=True)

        # Get FOV and Image size
        size, fov = BTAnalysisPipeline.get_osm_sky_dimensions("BTA/"+osm_file)

        print("Setting up OSKAR for "+osm_file)
        
        # FIXME: Implement config parser

        # Setup the interferometer ini file
        # FIXME: Write to new interferometer file and copy over settings
        #with open('sim_intif.ini', 'w', encoding='utf8'):
        
            # Set sky model location
            #ofname = "oskar_sky_model/file="+osm_file

            # Set frequency bin
            #freq = osm.split("_")[-1][:-7]+"e6"
            #ofname = r"start_frequency_hz="+freq
            #BTAnalysisPipeline.find_replace_line("BTA/test_intif.ini", "freqset", ofname)

        # Setup the imager ini file
        # FIXME: Write to new imager file and copy over settings
        #with open('sim_img.ini', 'w', encoding='utf8'):
            # Set the pixel resolution
            #ofname = "size="+str(size) 
            #BTAnalysisPipeline.find_replace_line("BTA/test_img.ini", "sizeset", ofname)

            # Set the FOV
            #ofname = "fov_deg="+str(fov)
            #BTAnalysisPipeline.find_replace_line("BTA/test_img.ini", "fovset", ofname)

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
    def setup_BTA_dir(h5_file, imager_template_ini="./regrid/test_intif_inis/test_img_gen.ini", interferometer_template_ini="./regrid/test_intif_inis/test_intif_gen.ini", oskar_telescope_model="./oskar_run_stage/telescope_model_AAstar", template=False):
        """
        Sets up the operating directory from which all anaysis will be done.

        :param h5_file: The location of the H5 file.
        :param cd_in: Whether to cd into the directory once finished or not.
        :param imager_template_ini: The file location of the OSKAR imager settings file.
        :param interferometer_template_ini: The file location of the OSKAR interferometer settings template file. The ini file must contain a line under [observation] with text "fset" in liu of the start_frequency_hz setting, and a line under [sky] with text "preset" in liu of the oskar_sky_model/file setting.
        :param oskar_telescope_model: The telescope model for OSKAR to use.
        :param template: If true, handle and return no h5 data.
        :return: The new location of the H5, imager ini, and interterometer template ini files, as well as the telescope model location. Also returns the PWD if cd_in is True.
        :param osm_dir: If a set of osms have already been generated, use the path to the directory specified to find the osms.
        """

        cwd = os.getcwd()

        # 1. Create Directory
        subprocess.run(["mkdir","-p","BTA"], check=True)

        # 2. Move H5 file and INIs to directory
        if not template: subprocess.run(["cp",h5_file,"BTA/analysis.h5"], check=True)
        subprocess.run(["cp",imager_template_ini,"BTA/imager_template.ini"], check=True)
        subprocess.run(["cp",interferometer_template_ini,"BTA/interferometer_template.ini"], check=True)
        if not os.path.isdir("BTA/telescope_model"):
            subprocess.run(["cp","-r",oskar_telescope_model,"BTA/telescope_model"], check=True)

        return ("BTA/analysis.h5" if not template else ""), "BTA/imager_template.ini", "BTA/interferometer_template.ini", "BTA/telescope_model", cwd

    @staticmethod
    def clean_BTA_dir(outdir, fits_output, clean=True):
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
    def H5_box_to_datacube(file, phase_ref_point = Regrid_Helper.ZENITH_530, require_regrid = True, max_freq_res = 100e6, imager_template_ini = "./regrid/test_intif_inis/test_img_gen.ini", interferometer_template_ini = "./regrid/test_intif_inis/test_intif_gen.ini", outdir = ".", clean = True, oskar_exec = Regrid_Helper.OSKAR_SIF, oskar_mode="singularity", oskar_telescope_model = Regrid_Helper.TELESCOPE, template_preset = "", coeval = True, load_osm=False):
        """
        Full pipeline function for transforming a H5 simulation box output into a FITS datacube.

        :param file: Location of the H5 file.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always regrid frequency bins, if false, regrid only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution in Hz.
        :param imager_template_ini: The file location of the OSKAR imager settings file.
        :param interferometer_template_ini: The file location of the OSKAR interferometer settings template file. The ini file must contain a line under [observation] with text "fset" in liu of the start_frequency_hz setting, and a line under [sky] with text "preset" in liu of the oskar_sky_model/file setting.
        :param outdir: The output location of the final FITS file.
        :param oskar_exec: The SIF file or location of compiled OSKAR binaries.
        :param oskar_mode: How shall OSKAR be run? Options include: python, binary, command, singularity.
        :param oskar_telescope_model: The telescope model for OSKAR to use.
        :param template_preset: Use a mock values array instead of a h5 file. Ignores any provided h5 file.
        :param coeval: If the H5 box is coeval or lightcone based.
        :param load_osm: If true load treat the file variable as if it were an OSM file.
        """

        # Expand paths
        if oskar_mode == "binary" or oskar_mode == "singularity":
            oskar_exec = Regrid_Helper.expand_path(oskar_exec)
        
        file = Regrid_Helper.expand_path(file)
        interferometer_template_ini = Regrid_Helper.expand_path(interferometer_template_ini)
        imager_template_ini = Regrid_Helper.expand_path(imager_template_ini)
        oskar_telescope_model = Regrid_Helper.expand_path(oskar_telescope_model)
        outdir = Regrid_Helper.expand_path(outdir)
        
        # Set templates
        if load_osm: template_preset = file.split('/')[-1][:-4]

        template_flag = (template_preset != "")

        print("Setting up BTA directory ...")
        h5_file, img_temp_ini, intif_temp_ini, _, _ = BTAnalysisPipeline.setup_BTA_dir(file, oskar_telescope_model=oskar_telescope_model, imager_template_ini=imager_template_ini, interferometer_template_ini=interferometer_template_ini, template=template_flag)

        # Create output file locations
        h5_id = file.split('/')[-1][:-3] if not template_flag else template_preset
        osm_output = Regrid_Helper.expand_path("BTA/" + h5_id + "_sky.osm")
        fits_output = Regrid_Helper.expand_path("BTA/" + h5_id + "_image.fits")
        
        # Run the OSM generator
        if not os.path.isfile(osm_output):
            if not load_osm:
                if template_flag:
                    print("Generating OSM files from template ...")
                    template_value = Regrid.mock_values(template_preset)
                    Regrid.generate_osm_from_simulation(template_value, osm_output=osm_output)
                else:
                    print("Generating OSM files from H5 ...")
                    Regrid.generate_osm_from_H5(h5_file, phase_ref_point=phase_ref_point, require_regrid=require_regrid, max_freq_res=max_freq_res, osm_output=osm_output, coeval=coeval)
            else:
                subprocess.run(["cp", file, osm_output], check=True)

        print("Running OSKAR ...")
        BTAnalysisPipeline.run_oskar_on_osms(osm_output, imager_template_ini=img_temp_ini, interferometer_template_ini=intif_temp_ini, fits_output=fits_output, oskar_exec=oskar_exec, oskar_mode=oskar_mode)

        # If clean is true remove all data relating to execution
        print("Cleaning up ...")
        BTAnalysisPipeline.clean_BTA_dir(outdir=outdir, fits_output=fits_output, clean=clean)


# Testing stage

#Regrid.generate_osm_from_H5("./regrid/yuxiang_bts/yuxiang1.h5", osm_output="./regrid/osm_output/yuxiang1_non_uniform_zenith_osm", coeval=True, require_regrid=False)
#Regrid.generate_osm_from_H5("./regrid/yuxiang_bts/yuxiang1.h5", osm_output="./regrid/osm_output/yuxiang1_00_osm", coeval=True, phase_ref_point=Regrid_Helper.ZERO_RADEC)
#Regrid.generate_osm_from_H5("./regrid/yuxiang_bts/yuxiang1.h5", osm_output="./regrid/osm_output/yuxiang1_zenith_osm", coeval=True)

#for template_preset in ["gaussian", "point", "random", "flat", "sinusoid", "point"]:
#    template_value = Regrid.mock_values(template_preset, scale=20)
#    Regrid.generate_osm_from_simulation(template_value, osm_output="./regrid/osm_output/"+template_preset+"_zenith_osm")

#Collator.collate_fits("./regrid/test_output/yuxiang1_fits", "./regrid/test_output")
#Collator.collate_fits("./regrid/test_output/yuxiangbad_fits", "./regrid/test_output")

#BTAnalysisPipeline.H5_box_to_datacube(None, template_preset="gaussian")

BTAnalysisPipeline.H5_box_to_datacube("./regrid/osm_output/yuxiang1_zenith_osm", oskar_exec=Regrid_Helper.OSKAR_BIN, load_osm=True, oskar_mode="binary", oskar_telescope_model=Regrid_Helper.TELESCOPE)