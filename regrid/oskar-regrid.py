#!/bin/python3
# pylint: disable=E1101
# pylint: disable=W0102

import numpy as np
import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM as fmodel, z_at_value as getz
from scipy.interpolate import make_interp_spline as misp
from astropy.coordinates import SkyCoord
#import copy
import os
import h5py
import pandas as pd

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
    def mock_values(preset, scale = 3, d = (100, 100, 100)):
        """
        Create an array of mock simulation values.

        :param preset: Mock brightness temperature array format. Options are {"flat", "random", "gaussian", "sinusoid"}
        :param scale: Define the Kelvin scale of the array (e.g. default value will create a uniform array of 3 K or a random array of 0 - 3 K).
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
        box_len = file.get('user_params').attrs['BOX_LEN'] * file.get('cosmo_params').attrs['hlittle']
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
    def generate_osm_from_simulation(values, voxels = None, d = (100, 100, 100), z_ref = 7, phase_ref_point = SkyCoord(ra=0*u.rad, dec=0*u.rad, frame='icrs'), require_regrid = True, max_freq_res = 100e6, v = (1, 1, 1), output_master_osm=False, osm_output="osm_output", cosmology=Cosmo()):
        """
        Generate a set of .osm files for an OSKAR sky model based on a Mpc**3 simulation output.

        :param values: The simulation datacube, must have shape d.
        :param voxels: An array describing a series of voxel dimensions corresponding to each simulation datacube element.
        :param d: Number of voxels in simulation in dimensions (x, y, t).
        :param z_ref: Refrence redshift, the ending redshift of the simulation.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always regrid frequency bins, if false, regrid only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution in Hz.
        :param v: If uniform spaxels = True, provides the initial voxel dimensions in h^-1 Mpc in dimensions (x, y, t), and auto-generates the voxel configuration array.
        :param output_master_osm: If true, output the entire datacube to one .osm file. If false, output a set of .osm files corresponding to each refrence frequency.
        :param osm_output: The directory to output the osm file(s) if output_master_osm is false.
        :param cosmology: The specific cosmology parameters in the form of a custom Cosmo object.
        """
        print("Initialising ...")

        # Calculate refrence comoving dist.
        Dz_ref = cosmology.z_to_Dz(z_ref)
        f_ref = cosmology.z_to_f(z_ref)

        # Set voxels array
        print("Creating mock voxels ...")
        voxels = np.full((*d, 3), v, dtype=np.float64)

        # Set regrid flag
        regrid_flag = require_regrid

        # Set cosmological redshift parameters
        Dz = Dz_ref
        Dz_val = Dz_ref.to_value(u.Mpc)
        z_prev = z_ref
        fq = f_ref.to_value(u.Hz)

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
                    dθ = dx/Dz_pix
                    dφ = dy/Dz_pix

                    Dz_val = Dz_val + dt # Increment the value of Dz by voxel dimension

                    # STEPS 2 & 3 - Convert line-of-sight comoving distance to frequency

                    df = 0

                    if x == 0 and y == 0:
                        # Determine corresponding redshifts
                        z_bot = z_prev
                        z_top = cosmology.Dz_to_z(Dz+dt*u.Mpc)

                        # Convert to frequency
                        f_bot = cosmology.z_to_f(z_bot)
                        f_top = cosmology.z_to_f(z_top)

                        # Store altered frequency bandwidth
                        df = np.abs((f_bot-f_top).to_value(u.Hz))

                        Dz = Dz + dt*u.Mpc # Increment the value of Dz by voxel dimension
                        z_prev = z_top # Top z in current box = bottom z in next box
                    else:
                        df = voxels[0, 0, t, 2]

                    # STEP 4 - Convert brightness temperature to pixel flux
                    fxy = fq + df/2
                      
                    # Calculate pixelated luminosity
                    Fv = Regrid.brightness_temperature_to_flux(Tb=Tb, fxy=fxy, dθ=dθ, dφ=dφ)

                    # Increment cumulative frequency
                    fq = fq + df

                    # STEPS 5 & 6 - Convert flat angular resolution to RA, Dec deviation
                    # Convert to l, m coordinates
                    dl = dθ / (2 * np.pi)
                    dm = dφ / (2 * np.pi)

                    # Convert l, m coordinates to spherical coords
                    dRA = np.arcsin(dl)
                    dDc = np.arcsin(dm/np.cos(np.arcsin(dm))) # FIXME Contact Jack Line over correctness of formula

                    # STEP 7.1 - Check if regridding is needed
                    if df > max_freq_res:
                        regrid_flag = regrid_flag or True

                    # Save values
                    voxels[x, y, t, 0] = dRA
                    voxels[x, y, t, 1] = dDc
                    voxels[x, y, t, 2] = df
                    values[x, y, t] = Fv

            print("\rTime step #", t, end="")

        print("\nTransforming complete.")

        # STEP 7 - Regrid frequency-dimension data if needed
        # FIXME Fix spline boundary conditions
        if regrid_flag:
            print("Performing regrid ...")

            for x in range(d[0]):
                for y in range(d[1]):

                    print("\rSpaxel # (", x, ",", y, ")", end="")

                    # Create interpolation B-spline
                    bspline = misp(np.cumsum(voxels[x, y, :, 2]), values[x, y, :], bc_type="natural")

                    # Perform regridding
                    new_freq, freq_bandw = np.linspace(0, np.sum(voxels[x, y, :, 2]), d[2], retstep=True)
                    new_flux = np.clip(bspline(new_freq), 0, None)

                    # Create array of uniform bin sizes
                    freq_bins = np.ones(d[2]) * freq_bandw

                    # Save variables
                    voxels[x, y, :, 2] = freq_bins
                    values[x, y, :] = new_flux

            print("\nRegrid complete.")

        else:
            print("No regrid required!")
        
        # STEP 8 - Write data to OSM file
        print("Configuring datacube for OSKAR file format ...")

        # RA, Dec centering function
        centering = (lambda x: (x - np.max(x)/2))

        # Cumulative sums are more important than voxel bins now, add half the mean to approximate centre
        rasum = centering(np.cumsum(voxels[:,:,:,0], axis=0)) + np.mean(voxels[:,:,:,0])/2
        decsum = centering(np.cumsum(voxels[:,:,:,1], axis=1)) + np.mean(voxels[:,:,:,1])/2
        freqsum = f_ref.to_value(u.Hz) - (np.cumsum(voxels[:,:,:,2], axis=2) + np.mean(voxels[:,:,:,2])/2)

        # Use Skycoords to calculate spherical RA, Dec offsets
        source_pos = phase_ref_point.spherical_offsets_by(rasum * u.rad, decsum * u.rad)
        RAs = source_pos.ra.to_value(u.deg)
        Dcs = source_pos.dec.to_value(u.deg)

        # Record data to file
        if output_master_osm:
            print("Recording data to .osm file")
            with open('reformatted.osm', 'w', encoding='utf8') as osm:

                # Clear file contents
                osm.truncate(0)

                # Add header lines
                osm.write("# Entries Key:\n")
                osm.write("#00.000000 +00.000000 0.0000+e00 0.0 0.0 0.0 000.000e6\n")
                osm.write("# RA       Dec        Stokes I   Q   U   V   Freq0\n")

                # Write OSM lines
                for x in range(d[0]):
                    for y in range(d[1]):
                        for t in range(d[2]):

                            # Format data
                            RAscn = np.char.zfill(np.format_float_positional(RAs[x, y, t], 6, False), 10)
                            Decln = np.char.zfill(np.format_float_positional(np.abs(Dcs[x, y, t]), 6, False), 9)
                            value = np.format_float_scientific(values[x, y, t], 4, False)
                            freq0 = np.format_float_positional(freqsum[x, y, t] / 1e6, 3, False)

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
                                "\n"
                            )

                        print("\rSpaxel # (", x, ",", y, ")", end="")

            print("\nProcess complete, data saved to "+osm_output)
        else:
            print("Recording data to .osm files")

            if not os.path.isdir(osm_output):
                os.mkdir(osm_output)

            for t in range(d[2]):
                file_freq = np.format_float_positional(freqsum[0, 0, t] / 1e6, 3, False)
                print("\rGenerating OSM for freq0 =", file_freq, "MHz ( file #", t+1, "of", d[2], ")", end="")

                with open(osm_output+'/reformatted_no.'+str(t+1)+'_'+str(file_freq)+'MHz.osm', 'w', encoding='utf8') as osm:
                    # Clear file contents
                    osm.truncate(0)

                    # Add header lines
                    osm.write("# Entries Key:\n")
                    osm.write("#00.000000 +00.000000 0.0000+e00 0.0 0.0 0.0 000.000e6\n")
                    osm.write("# RA       Dec        Stokes I   Q   U   V   Freq0\n")

                    # Write OSM lines
                    for x in range(d[0]):
                        for y in range(d[1]):
                            # Format data
                            RAscn = np.char.zfill(np.format_float_positional(RAs[x, y, t], 6, False), 10)
                            Decln = np.char.zfill(np.format_float_positional(np.abs(Dcs[x, y, t]), 6, False), 9)
                            value = np.format_float_scientific(values[x, y, t], 4, False)
                            freq0 = np.format_float_positional(freqsum[x, y, t] / 1e6, 3, False)

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
                                "\n"
                            )

            print("\nProcess complete, data saved to "+osm_output)
    
    @staticmethod
    def generate_osm_from_H5(file, phase_ref_point = SkyCoord(ra=0*u.rad, dec=0*u.rad, frame='icrs'), require_regrid = True, max_freq_res = 100e6, output_master_osm=False, osm_output="", coeval=True):
        """
        Combines both the convert_H5_to_csv and generate_osm_from_simulation functions.

        :param file: Location of the H5 file.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always regrid frequency bins, if false, regrid only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution in Hz.
        :param output_master_osm: If true, output the entire datacube to one .osm file. If false, output a set of .osm files corresponding to each refrence frequency.
        :param osm_output: The directory to output the osm file(s) if output_master_osm is false.
        """

        values, dim, z_ref, vox, cosmology = Regrid.convert_H5_to_csv(file, coeval=coeval)

        if osm_output == "": osm_output = file.split('/')[-1][:-3] + "_osm"

        Regrid.generate_osm_from_simulation(values, d=dim, z_ref=z_ref, require_regrid=require_regrid, max_freq_res=max_freq_res, v=vox, output_master_osm=output_master_osm, osm_output=osm_output, cosmology=cosmology, phase_ref_point=phase_ref_point)

class Collator(object):
    """
    The Collator class contains methods pertaining to the collation of output FITS file images into a spectral datacube.
    """

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
    def collate_header(h5_file="", osm_dir="", fits_dir=""):
        """
        Compiles a multiline string containing a functional FITS header for a collated datacube. Any parameter set to none will not have information taken from it.

        :param h5_file: The original simulation H5 file.
        :param osm_dir: The outputted OSM file directory.
        :param fits_dir: The outputted FITS file directory.
        :return: The file header as a dictionary.
        """

        # FIXME Add a header modification for outname.

        header = {}

        if h5_file == "":
            print("Retreiving H5 header information ...")
        if osm_dir == "":
            print("Retreiving OSM header information ...")
        if fits_dir == "":
            print("Retreiving FITS header information ...")

        return header

    @staticmethod
    def collate_fits(fits_dir, outdir=".", headers=["", "", ""]):
        """
        Collates a directory of FITS images into a single FITS datacube. Adds a header to the file that corresponds to useful information about the simulation box.

        :param fits_dir: Directory containing the FITS files to be collated.
        :param outdir: The output directory of the FITS file.
        :param headers: A list containing the H5 file, the OSM directory, and FITS directory locations for header compilation.
        :return: The file name of the FITS datacube.
        """
        print("Collating")

        
        img_list = []

        fits_files = Collator.dir_list_sorted(dir_=fits_dir)

        for file in fits_files:
            img_list.append(fits.getdata(fits_dir+"/"+file))
            print("\rCollating file # (", file, ")", end="")
        
        img_arr = np.array(img_list)
        outname = outdir+"/"+fits_dir.split("/")[-1].split(".")[0] + "_cube.fits"

        hdu_new = fits.PrimaryHDU(img_arr)

        print("Retreiving header information ...")
        header = Collator.collate_header(*headers)

        for item in header.keys():
            hdu_new.header[item] = header[item]

        print("\nSaving to file: " + outname)
        hdu_new.writeto(outname)

        print("Collation Complete")

        return outname

class BTAnalysisPipeline(object):
    """
    A broader class that combines all components of the individual components of the simulated IGM to simulated observation pipeline together.
    """

    @staticmethod
    def get_osm_sky_dpi(osm_file, fov=5):
        """
        Get the pixel count for a FITS file with it's corresponding OSKAR sky model.

        :param osm_file: The OSM file to analyse.
        :return: The amount of pixels in per FOV.
        """

        df = pd.read_csv(osm_file, delimiter=" ", skiprows=3, index_col=False, names=["RA", "Dec", "Stokes I", "Q", "U", "V", "Freq0"])

        # Calculate RA dimension
        RAc = np.abs(np.array(df['RA'])[-1] - df['RA'][0]) % 360

        # Calculate Dec dimension
        Decc = np.abs(np.array(df['Dec'])[-1] - df['Dec'][0] + 360) % 360

        # Get number of voxels on a side
        n = np.sqrt(len(np.array(df['RA'])))

        # Calculate pixel length of FOV image
        dpi = n * fov / max(RAc, Decc)

        return dpi

    @staticmethod
    def run_oskar_on_osms(osm_dir, imager_template_ini = "./regrid/test_intif_inis/test_img_gen.ini", interferometer_template_ini = "./test_intif_inis/test_intif_gen.ini", fits_output="", oskar_sif="~/.oskar/OSKAR-2.8.3-Python3.sif", fov=5.0):
        """
        Run oskar on each of the OSM sky models found in a fits directory, should already be formatted according to the output of the Regrid object.

        :param osm_dir: Directory containing the OSM files to be imaged.
        :param imager_template_ini: The file location of the OSKAR imager settings file.
        :param interferometer_template_ini: The file location of the OSKAR interferometer settings template file. The ini file must contain a line under [observation] with text "fset" in liu of the start_frequency_hz setting, and a line under [sky] with text "preset" in liu of the oskar_sky_model/file setting.
        :param fits_output: The directory to output the resultant FITS files.
        :param oskar_sif: The SIF file containing the OSKAR program. OSKAR must be run from singularity.
        """
        # FIXME: Run OSKAR Python Interface

        # Create the output file
        os.system("mkdir -p BTA/output")

        # Get files to iterate over
        osm_list = Collator.dir_list_sorted(dir_="BTA/"+osm_dir)

        for osm in osm_list:
            # Setup the interferometer ini file
            # Duplicate the generator file
            os.system("cp "+interferometer_template_ini+" BTA/test_intif.ini")

            # Set sky model location
            ofname = r"oskar_sky_model\/file=BTA\/"+osm_dir+r"\/"+osm
            os.system(r'sed -i "s/^preset.*/'+ofname+r'/" "BTA/test_intif.ini"')

            # Set frequency bin
            freq = osm.split("_")[-1][:-7]+"e6"
            ofname = r"start_frequency_hz="+freq
            os.system(r'sed -i "s/^fset.*/'+ofname+r'/" "BTA/test_intif.ini"')

            # Run OSKAR's interferometer simulation
            os.system("singularity exec --nv --bind $PWD --cleanenv --home $PWD "+oskar_sif+" oskar_sim_interferometer BTA/test_intif.ini")

            # Setup the imager ini file
            # Duplicate the generator file
            os.system("cp "+imager_template_ini+" BTA/test_img.ini")

            # Set the pixel resolution
            size = BTAnalysisPipeline.get_osm_sky_dpi("BTA/"+osm_dir+"/"+osm, fov=fov)
            ofname = "size="+str(size) 
            os.system(r'sed -i "s/^sizeset.*/'+ofname+r'/" "BTA/test_img.ini"')

            # Set the FOV
            ofname = "fov_deg="+str(fov)
            os.system(r'sed -i "s/^fovset.*/'+ofname+r'/" "BTA/test_img.ini"')

            # Run OSKAR's imager simulation
            os.system("singularity exec --nv --bind $PWD --cleanenv --home $PWD "+oskar_sif+" oskar_sim_interferometer BTA/test_intif.ini")


    @staticmethod
    def setup_BTA_dir(h5_file, cd_in=True, imager_template_ini="./regrid/test_intif_inis/test_img_gen.ini", interferometer_template_ini="./regrid/test_intif_inis/test_intif_gen.ini", oskar_telescope_model="~/.oskar/telescope_model_AAstar", template=False):
        """
        Sets up the operating directory from which all anaysis will be done.

        :param h5_file: The location of the H5 file.
        :param cd_in: Whether to cd into the directory once finished or not.
        :param imager_template_ini: The file location of the OSKAR imager settings file.
        :param interferometer_template_ini: The file location of the OSKAR interferometer settings template file. The ini file must contain a line under [observation] with text "fset" in liu of the start_frequency_hz setting, and a line under [sky] with text "preset" in liu of the oskar_sky_model/file setting.
        :param oskar_telescope_model: The telescope model for OSKAR to use.
        :param template: If true, handle and return no h5 data.
        :return: The new location of the H5, imager ini, and interterometer template ini files, as well as the telescope model location. Also returns the PWD if cd_in is True.
        """

        cwd = os.getcwd()

        # 1. create directory
        os.system("mkdir -p BTA")

        # 2. Move H5 file and INIs to directory
        if not template: os.system("cp "+h5_file+" BTA/analysis.h5")
        os.system("cp "+imager_template_ini+" BTA/imager_template.ini")
        os.system("cp "+interferometer_template_ini+" BTA/interferometer_template.ini")
        os.system("cp -r "+oskar_telescope_model+" BTA/telescope_model")

        # 3. CD into directory
        if cd_in:
            os.system("cd BTA")
            return ("BTA/analysis.h5" if not template else ""), "BTA/imager_template.ini", "BTA/interferometer_template.ini", "BTA/telescope_model", cwd

    @staticmethod
    def clean_BTA_dir(outdir, fits_cube, cd_out=True, clean=True):
        """
        Finishes and cleans up the mess created by the BTA class.

        :param outdir: The directory to save the finalised fits image, relative to parent dir.
        :param fits_cube: The name of the FITS output file.
        :param cd_out: Whether to cd out to parent directory once done.
        :param clean: Whether or not to remove the BTA directory. Only works if cd_out is true.
        """

        # Move datacube to parent directory
        os.system("mv "+fits_cube+" ../"+outdir+"/"+fits_cube)

        # CD to parent dir
        if cd_out:
            os.system("cd ..")

            if clean:
                os.system("rm -rf BTA")

    @staticmethod
    def H5_box_to_datacube(file, phase_ref_point = SkyCoord(ra=0*u.rad, dec=0*u.rad, frame='icrs'), require_regrid = True, max_freq_res = 100e6, imager_template_ini = "./regrid/test_intif_inis/test_img_gen.ini", interferometer_template_ini = "./regrid/test_intif_inis/test_intif_gen.ini", outdir = ".", clean=True, oskar_sif="~/.oskar/OSKAR-2.8.3-Python3.sif", oskar_telescope_model="~/.oskar/telescope_model_AAstar", fov=5.0, template_preset="", coeval=True):
        """
        Full pipeline function for transforming a H5 simulation box output into a FITS datacube.

        :param file: Location of the H5 file.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always regrid frequency bins, if false, regrid only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution in Hz.
        :param imager_template_ini: The file location of the OSKAR imager settings file.
        :param interferometer_template_ini: The file location of the OSKAR interferometer settings template file. The ini file must contain a line under [observation] with text "fset" in liu of the start_frequency_hz setting, and a line under [sky] with text "preset" in liu of the oskar_sky_model/file setting.
        :param outdir: The output location of the final FITS file.
        :param oskar_sif: The SIF file containing the OSKAR program. OSKAR must be run from singularity.
        :param oskar_telescope_model: The telescope model for OSKAR to use.
        :param fov: The field of view of imaging/observation
        :param template_preset: Use a mock values array instead of a h5 file. Ignores any provided h5 file.
        """

        template_flag = (template_preset != "")

        print("Setting up BTA directory ...")
        h5_file, img_temp_ini, intif_temp_ini, _, _ = BTAnalysisPipeline.setup_BTA_dir(file, oskar_telescope_model=oskar_telescope_model, imager_template_ini=imager_template_ini, interferometer_template_ini=interferometer_template_ini, template=template_flag)

        osm_output = (file.split('/')[-1][:-3] if not template_flag else template_preset) + "_osm"
        fits_output = (file.split('/')[-1][:-3] if not template_flag else template_preset) + "_fits"
        
        if template_flag:
            print("Generating OSM files from template ...")
            template_value = Regrid.mock_values(template_preset)
            Regrid.generate_osm_from_simulation(template_value, osm_output="BTA/"+osm_output)
        else:
            print("Generating OSM files from H5 ...")
            Regrid.generate_osm_from_H5(h5_file, phase_ref_point=phase_ref_point, require_regrid=require_regrid, max_freq_res=max_freq_res, osm_output="BTA/"+osm_output, coeval=coeval)

        print("Running OSKAR. Outputting to ./BTA/oskar.out ...")
        BTAnalysisPipeline.run_oskar_on_osms(osm_output, imager_template_ini=img_temp_ini, interferometer_template_ini=intif_temp_ini, fits_output=fits_output, oskar_sif=oskar_sif, fov=fov)

        print("Collating fits images ...")
        fits_cube = Collator.collate_fits(fits_output, headers=[h5_file, osm_output, fits_output])

        # If clean is true remove all data relating to execution
        print("Cleaning up ...")
        BTAnalysisPipeline.clean_BTA_dir(outdir=outdir, fits_cube=fits_cube, clean=clean)


# Testing stage

Regrid.generate_osm_from_H5("./regrid/yuxiang_bts/yuxiang1.h5", osm_output="./regrid/yuxiang1_osm", coeval=True)

#Collator.collate_fits("./regrid/test_output/yuxiang1_fits", "./regrid/test_output")

#BTAnalysisPipeline.H5_box_to_datacube(None, template_preset="gaussian")