#!/bin/python3
# pylint: disable=E1101

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
    def convert_H5_to_csv(h5_location, save_data=False, outdir='', name="out_h5_data"):
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
    def generate_osm_from_simulation(values, voxels = None, d = (100, 100, 100), z_ref = 7, phase_ref_point = SkyCoord(ra=0*u.rad, dec=0*u.rad, frame='icrs'), require_regrid = True, max_freq_res = 100e6, uniform_spaxels = True, v = (1, 1, 1), output_master_osm=False, osm_output="osm_output", cosmology=Cosmo()):
        """
        Generate a set of .osm files for an OSKAR sky model based on a Mpc**3 simulation output.

        :param values: The simulation datacube, must have shape d.
        :param voxels: An array describing a series of voxel dimensions corresponding to each simulation datacube element.
        :param d: Number of voxels in simulation in dimensions (x, y, t).
        :param z_ref: Refrence redshift, the ending redshift of the simulation.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always regrid frequency bins, if false, regrid only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution in Hz.
        :param uniform_spaxels: If true then each voxel has the same physical dimensions, the dimensions are given by v. The voxel array is automatically generated.
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
        if uniform_spaxels:
            print("Creating mock voxels ...")
            
            voxels = np.full((*d, 3), v, dtype=np.float64)

        # Set regrid flag
        regrid_flag = require_regrid

        print("Transforming coordinates ...")
        # Main loop of creation
        for t in range(d[2]):
            for x in range(d[0]):
                for y in range(d[1]):

                    Dz = Dz_ref
                    Dz_val = Dz_ref.to_value(u.Mpc)
                    z_prev = z_ref
                    fq = f_ref.to_value(u.Hz)

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

                    if (not uniform_spaxels) or (x == 0 and y == 0):
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
                    Fv = 2 * c.k_B.value * fxy**2 * Tb * (dθ * dφ) / c.c.value ** 2
                    Fv = Fv * 1e26 # Converts to Jansky

                    # Increment cumulative frequency
                    fq = fq + df

                    # STEPS 5 & 6 - Convert flat angular resolution to RA, Dec deviation
                    # Convert to l, m coordinates
                    dl = dθ / (2 * np.pi)
                    dm = dφ / (2 * np.pi)

                    # Convert l, m coordinates to spherical coords
                    dRA = np.arcsin(dl)
                    dDc = np.arcsin(dm/np.cos(np.arcsin(dm)))

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
        if regrid_flag:
            print("Performing regrid ...")

            for x in range(d[0]):
                for y in range(d[1]):

                    print("\rSpaxel # (", x, ",", y, ")", end="")

                    # Create interpolation B-spline
                    bspline = misp(np.cumsum(voxels[x, y, :, 2]), values[x, y, :])

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
        rasum = centering(np.cumsum(voxels[:,:,:,0], axis=0)) + np.mean(voxels[:,:,:,0], axis=0)/2
        decsum = centering(np.cumsum(voxels[:,:,:,1], axis=1)) + np.mean(voxels[:,:,:,1], axis=1)/2
        freqsum = f_ref.to_value(u.Hz) - np.cumsum(voxels[:,:,:,2], axis=2) - np.mean(voxels[:,:,:,2], axis=2)/2

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

            print("\nProcess complete, data saved to ./reformatted.osm")
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

            print("\nProcess complete, data saved to ./osm-output/")
    
    @staticmethod
    def generate_osm_from_H5(file, phase_ref_point = SkyCoord(ra=0*u.rad, dec=0*u.rad, frame='icrs'), require_regrid = True, max_freq_res = 100e6, uniform_spaxels = True, output_master_osm=False):
        """
        Combines both the convert_H5_to_csv and generate_osm_from_simulation functions.

        :param file: Location of the H5 file.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always regrid frequency bins, if false, regrid only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution in Hz.
        :param uniform_spaxels: If true then each voxel has the same physical dimensions, the dimensions are given by v. The voxel array is automatically generated.
        :param output_master_osm: If true, output the entire datacube to one .osm file. If false, output a set of .osm files corresponding to each refrence frequency.
        """

        values, dim, z_ref, vox, cosmology = Regrid.convert_H5_to_csv(file)

        osm_output = file.split('/')[-1][:-3] + "_osm"

        Regrid.generate_osm_from_simulation(values, d=dim, z_ref=z_ref, require_regrid=require_regrid, max_freq_res=max_freq_res, uniform_spaxels=uniform_spaxels, v=vox, output_master_osm=output_master_osm, osm_output=osm_output, cosmology=cosmology, phase_ref_point=phase_ref_point)

class Collator(object):
    """
    The Collator class contains methods pertaining to the collation of output FITS file images into a spectral datacube.
    """

    @staticmethod
    def collate_fits(fits_dir, outdir="."):
        """
        Collates a directory of FITS images into a single FITS datacube. Adds a header to the file that corresponds to useful information about the simulation box.

        :param fits_dir: Directory containing the FITS files to be collated.
        :return: The file name of the FITS datacube.
        """
        print("Collating")

        fits_files = os.listdir(fits_dir)
        img_list = []

        for file in fits_files:
            img_list.append(fits.getdata(fits_dir+"/"+file))
            print("\rCollating file # (", file, ")", end="")
        
        img_arr = np.array(img_list)
        outname = outdir+"/"+fits_dir.split("/")[-1].split(".")[0] + "_cube.fits"

        print("\nSaving to file: " + outname)

        hdu_new = fits.PrimaryHDU(img_arr)
        hdu_new.writeto(outname)

        print("Collation Complete")

        # FIXME Add a header modification for outname.

        return outname

class BTAnalysisPipeline(object):
    """
    A broader class that combines all components of the individual components of the simulated IGM to simulated observation pipeline together.
    """

    @staticmethod
    def run_oskar_on_osms(osm_dir, imager_ini = "../oskar_run_stage/oskar_imager.ini", interferometer_template_ini = "./test_intif_inis/test_intif_gen.ini"):
        """
        Run oskar on each of the OSM sky models found in a fits directory, should already be formatted according to the output of the Regrid object.

        :param osm_dir: Directory containing the OSM files to be imaged.
        :param imager_ini: The file location of the OSKAR imager settings file.
        :param interferometer_template_ini: The file location of the OSKAR interferometer settings template file. The ini file must contain a line under [observation] with text "fset" in liu of the start_frequency_hz setting, and a line under [sky] with text "preset" in liu of the oskar_sky_model/file setting.
        """

    @staticmethod
    def setup_BTA_dir(h5_file, cd_in=True):
        """
        Sets up the operating directory from which all anaysis will be done.

        :param h5_file: The location of the H5 file.
        :param cd_in: Whether to cd into the directory once finished or not.
        :param imager_ini: The file location of the OSKAR imager settings file.
        :param interferometer_template_ini: The file location of the OSKAR interferometer settings template file. The ini file must contain a line under [observation] with text "fset" in liu of the start_frequency_hz setting, and a line under [sky] with text "preset" in liu of the oskar_sky_model/file setting.
        :return: The new location of the H5, imager ini, and interterometer template ini files. Also returns the PWD if cd_in is True.
        """

        # 1. create directory
        # 2. Move H5 file to directory
        # 3. CD into directory

        return "", "", "", ""

    @staticmethod
    def clean_BTA_dir(pwd, outdir, fits_cube, clean=True):
        """
        Finishes and cleans up the mess created by the BTA class.

        :param pwd: The directory to cd to once cleaning is complete.
        :param outdir: The directory to save the finalised fits image, relative to pwd.
        :param fits_cube: The name of the FITS output file.
        :param clean: Whether or not to remove the BTA directory. 
        """

    @staticmethod
    def H5_box_to_datacube(file, phase_ref_point = SkyCoord(ra=0*u.rad, dec=0*u.rad, frame='icrs'), require_regrid = True, max_freq_res = 100e6, uniform_spaxels = True, imager_ini = "../oskar_run_stage/oskar_imager.ini", interferometer_template_ini = "./test_intif_inis/test_intif_gen.ini", outdir = ".", clean=True):
        """
        Full pipeline function for transforming a H5 simulation box output into a FITS datacube.

        :param file: Location of the H5 file.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always regrid frequency bins, if false, regrid only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution in Hz.
        :param uniform_spaxels: If true then each voxel has the same physical dimensions, the dimensions are given by v. The voxel array is automatically generated.
        :param imager_ini: The file location of the OSKAR imager settings file.
        :param interferometer_template_ini: The file location of the OSKAR interferometer settings template file. The ini file must contain a line under [observation] with text "fset" in liu of the start_frequency_hz setting, and a line under [sky] with text "preset" in liu of the oskar_sky_model/file setting.
        :param outdir: The output location of the final FITS file.
        """

        print("Setting up BTA directory ...")
        h5_file, img_ini, intif_temp_ini, pwd = BTAnalysisPipeline.setup_BTA_dir(file)
        osm_output = file.split('/')[-1][:-3] + "_osm"
        fits_output = file.split('/')[-1][:-3] + "_fits"
        
        print("Generating OSM files from H5 ...")
        Regrid.generate_osm_from_H5(h5_file, phase_ref_point=phase_ref_point, require_regrid=require_regrid, max_freq_res=max_freq_res, uniform_spaxels=uniform_spaxels)

        print("Running OSKAR. Outputting to ./BTA/oskar.out ...")
        BTAnalysisPipeline.run_oskar_on_osms(osm_output, imager_ini=img_ini, interferometer_template_ini=intif_temp_ini)

        print("Collating fits images ...")
        fits_cube = Collator.collate_fits(fits_output)

        # If clean is true remove all data relating to execution
        print("Cleaning up ...")
        BTAnalysisPipeline.clean_BTA_dir(pwd=pwd, outdir=outdir, fits_cube=fits_cube, clean=clean)


# Testing stage

for preset_opt in ("point", "gaussian", "yuxiang1", "yuxiang2"):
    Collator.collate_fits("./regrid/test_output/"+preset_opt+"_fits", outdir="./regrid/test_output")