"""
The reformatter module contains functions relating to translating simulation data to OSKAR output data.
"""

# System imports
import warnings

# Mathematics and calculations
import astropy.constants as c
import astropy.units as u
import numpy as np
import h5py

# Data handling and statistics
import pandas as pd
from scipy.interpolate import make_interp_spline as misp

# Astropy extras
from astropy.coordinates import SkyCoord

# Local imports
from oskareor.skalow_calc import EoRCosmology as eorcosmo, SKAMath as omath, OSKARFileConfig as ofc
from oskareor.oskar_helpers import OSKARHelper as ohelp

class SimulationReformatter():
    """
    The reformatter class contains functions relating to translating simulation data to OSKAR output data.

    NB: all units are assumed SI unless explicitly stated as an astropy quantity. Exceptions to this rule include:
    - Flux (in Jansky)
    - Distance (in Mpc)
    
    Sky coordinate locations are given in Degrees, whereas angular diameters are given in radians.
    """

    @staticmethod
    def brightness_temperature_to_flux(tb: float, fxy: float, dtheta: float, dphi: float) -> float:
        """
        Given the temperature, frequency, and angular size in both sky directions, calculate the luminosity in Jy.
        
        :param tb: The brightness temperature. Must be provided in Kelvin.
        :param fxy: The frequency of emission. Must be provided in units of Hertz.
        :param dtheta: The first specified angular diameter. Perpendicular to the first. Must be provided in units of Radians.
        :param dphi: The second specified angular diameter. Perpendicular to the first. Must be provided in units of Radians.

        :return fv: The flux corresponding to the brightness temperature. In units of Jansky.
        """

        # Calculate pixelated luminosity
        fv = 2 * c.k_B.si.value * fxy**2 * tb * (dtheta * dphi) / c.c.si.value ** 2

        return fv * 1e26 # Converts to Jansky
    
    @staticmethod
    def brightness_temperature_to_linewidth(tb: float) -> float:
        """
        Takes the brightness temperature emission. Assumes that the kinetic gas temperature is coupled with the spin temperature and the brightness temperature (generally correct for z < 25).

        :param tb: The brightness temperature. Must be provided in Kelvin.

        :return linewidth: The flux corresponding to the brightness temperature. In units of Jansky.
        """

        # Calculate linewidth in Hz
        return (tb ** 0.5) * omath.SIGMA_F_FLOAT
    
    TEMPLATE_PRESETS = {
        "gaussian" : (
            { "normal", "gauss", "bell", "n", "g", "b" },
            "Rotationally symmetric centered gaussian plane with FWHM = d(t)/3.",
            lambda p: omath.gaussian(p['r'], var=omath.STDEV(p['d'][2]/3)**2, amp=p['T_max'])
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
            lambda p: np.abs(omath.sinc(p['r'], f=1/p['d'][2], amp=p['T_max']))
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
            lambda p: ohelp.COEVAL_TEMPLATE_VALUES_1[*((200 + ((np.array([p['i'], p['j'], p['t']]) - (np.array(p['d']) // 2)))) % 400)]
            ),
        "coeval2"  : (
            { "coeval 2", "2", "yuxiang2", "yuxiang 2", "y2", "c2" },
            "The second of two simulation boxes, cocentric with the desired values box, the original model has d = (400, 400, 400).\n"
            + "If d(a) < 400 then the box outer edges will be cropped and if d(a) > 400 the box will repeat beyond 400 px from the centre.",
            lambda p: ohelp.COEVAL_TEMPLATE_VALUES_2[*((200 + ((np.array([p['i'], p['j'], p['t']]) - (np.array(p['d']) // 2)))) % 400)]
            )
    }

    @staticmethod
    def display_template_presets(print_presets: bool = True, filter_preset: str = "") -> dict:
        """
        Return and/or display all available templates, their names, and their descriptions. All templates except random, coeval1, and coeval2 are identical in all t-dimension voxels.

        :param print_presets: If true print the templates to console and return the dictionary. If false only return the dictionary.
        :param filter_preset: If true print/return the information for only one specific template entry, according to the given string.
        
        :return template_options: The dictionary containing all available templates or a dictionary of the specific desired template.
        """

        return ohelp.display_options(SimulationReformatter.TEMPLATE_PRESETS, print_options=print_presets, selection=filter_preset)

    @staticmethod
    def mock_values(preset: str, scale: float = 10, d: tuple = (100, 100, 100), special = None):
        """
        Create an array of mock simulation values.

        :param preset: Mock brightness temperature array format. Run SimulationReformatter.display_template_presets for more information.
        :param scale: a.k.a. `T_max`. The maximum Kelvin value for the whole array, acts as a normalisation factor.
        :param d: The size of the values datacube.
        :param special: A custom lambda function that takes the dictionary of parameters (`d`, `i`, `j`, `x`, `y`, `t`, `r`, `T_max`) and returns a float, treat the preset parameter as a custom name.
        
        Note that `x` and `y` are positioned so that the centermost pixel is (0, 0) whereas `i` and `j` are the standard array values array indicies. Only `d`, `i`, and `j` are indicies, the others should be treated as floats.

        :return values: Mock brightness temperature values.
        """

        # Select specific template
        if special is None:
            selection = SimulationReformatter.display_template_presets(False, preset)
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

            print("\rCreating values for t-dim slice #", t, "of", d[2], end="")
        
        print("\nValues created!")

        return values.astype(np.float64)

    @staticmethod
    def convert_h5_coeval_to_csv(h5_location, save_data = False, outdir = '', name = "out_h5_data"):
        """
        Extract h5 data from Yuxiang Qin's simulations and either return the data objects or save to a seperate CSV.

        :param h5_location: The file location of the h5 data.
        :param save_data: If true, output data to a CSV and text file, specified by the outdir parameter.
        :param outdir: The directory to output both CSV and text information.
        :param name: The file name template to be saved to.

        :return (values, d, z_ref, vox, cosmo): The numpy values array in Kelvin, the shape of the array, the refrence redshift, the voxel size in Mpc, and the simulation box cosmology.
        """

        file = h5py.File(h5_location, 'r')

        # Get BT data
        bt_data = np.array(file.get('BrightnessTemp')['brightness_temp'])

        # Get Box and Voxel dimensions
        box_len = file.get('user_params').attrs['BOX_LEN']#* file.get('cosmo_params').attrs['hlittle']
        vox = np.array(np.ones(3) * box_len / bt_data.shape[0])

        # Define cosmology with H0=100h
        cosmology = eorcosmo(
            omega_m_0 = file.get('cosmo_params').attrs['OMm'],
            omega_b_0 = file.get('cosmo_params').attrs['OMb']
        )

        # Transform intitial redshift
        z_mid = file.attrs['redshift']
        z_ref = cosmology.dz_to_z(cosmology.z_to_dz(z_mid) - u.Mpc*box_len/2)

        if save_data:
            np.savetxt(outdir+'/'+name+'.csv', bt_data, delimiter=", ")
            np.savetxt(outdir+'/'+name+'.txt', np.array([bt_data.shape, z_ref, vox, cosmology]), delimiter=", ")

        return bt_data, bt_data.shape, z_ref, vox, cosmology
    
    @staticmethod
    def convert_h5_lightcone_to_csv(h5_location, save_data = False, outdir = '', name = "out_h5_data"):
        """
        Extract h5 data from Yuxiang Qin's lightcone simulations and either return the data objects or save to a seperate CSV.

        :param h5_location: The file location of the h5 data.
        :param save_data: If true, output data to a CSV and text file, specified by the outdir parameter.
        :param outdir: The directory to output both CSV and text information.
        :param name: The file name template to be saved to.
        :return (values, d, z_ref, vox, cosmo): The numpy values array in Kelvin, the shape of the array, the refrence redshift, the voxel size in Mpc, and the simulation box cosmology.
        """

        file = h5py.File(h5_location, 'r')

        # Get BT data
        bt_data = np.array(file.get('lightcones/brightness_temp'))

        # Define cosmology with H0=100h
        cosmology = eorcosmo(
            omega_m_0 = file.get('cosmo_params').attrs['OMm'],
            omega_b_0 = file.get('cosmo_params').attrs['OMb']
        )

        # Get Box and Voxel dimensions
        box_len = file.get('user_params').attrs['BOX_LEN'] * file.get('cosmo_params').attrs['hlittle']
        rs = np.array(list(file.get('node_redshifts')))
        vox = np.array([box_len, box_len, abs(cosmology.z_to_dz(rs[0])-cosmology.z_to_dz(rs[-1])).to_value(u.Mpc)]) / bt_data.shape
        
        # Transform intitial redshift
        z_ref = np.min(rs)

        if save_data:
            np.savetxt(outdir+'/'+name+'.csv', bt_data, delimiter=", ")
            np.savetxt(outdir+'/'+name+'.txt', np.array([bt_data.shape, z_ref, vox, cosmology]), delimiter=", ")

        return bt_data, bt_data.shape, z_ref, vox, cosmology
    
    @staticmethod
    def convert_h5_to_csv(h5_location: str, save_data: bool = False, outdir: str = '', name: str = "out_h5_data", coeval: bool = True):
        """
        Extract h5 data from Yuxiang Qin's simulations and either return the data objects or save to a seperate CSV.

        :param h5_location: The file location of the h5 data.
        :param save_data: If true, output data to a CSV and text file, specified by the outdir parameter.
        :param outdir: The directory to output both CSV and text information.
        :param name: The file name template to be saved to.
        :param coeval: whether or not the box is coeval or lightcone.
        :return (values, d, z_ref, vox, cosmo): The numpy values array in Kelvin, the shape of the array, the refrence redshift, the voxel size in Mpc, and the simulation box cosmology.
        """

        if coeval:
            return SimulationReformatter.convert_h5_coeval_to_csv(h5_location=h5_location, save_data=save_data, outdir=outdir, name=name)
        
        return SimulationReformatter.convert_h5_lightcone_to_csv(h5_location=h5_location, save_data=save_data, outdir=outdir, name=name)
        
    @staticmethod
    def transform_datacube_units(values, voxels, z_ref = 7, require_regrid = True, max_freq_res = 100 * u.MHz, cosmology = eorcosmo()):
        """
        Transform a datacube with dimensions x, y, t (cMpc x cMpc x cMpc) to ⍺, δ, f (rad x rad x Hz),

        :param values: The simulation datacube in units of Kelvin.
        :param voxels: An array describing a series of voxel dimensions corresponding to each simulation datacube voxel element in units of (cMpc, cMpc, cMps).
        :param z_ref: Refrence redshift, the ending redshift of the simulation.
        :param require_regrid: If true then always reformat frequency bins, if false, reformat only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution.
        :param v: If all voxels are the same, provides the initial voxel dimensions in cMpc in dimensions (x, y, t), and auto-generates the voxel configuration array.
        :param osm_output: The relative path to save the osm file to.
        :param cosmology: The specific cosmology parameters in the form of a custom eorcosmo object.
        """

        # Configure d variable
        d = np.shape(values)

        # Set regrid flag
        regrid_flag = require_regrid

        # Calculate refrence comoving dist.
        dz_ref = cosmology.z_to_dz(z_ref)
        f_ref = cosmology.z_to_f(z_ref)

        # Set cosmological redshift parameters
        dz = dz_ref
        z_prev = z_ref
        fq = f_ref.to_value(u.Hz)
        dz_val = dz_ref.to_value(u.Mpc)

        # Set maximum frequency resolution
        max_freq_res_hz = max_freq_res.to_value(u.Hz)

        # Set Linewidth array
        sigma_f = SimulationReformatter.brightness_temperature_to_linewidth(values)

        print("Transforming coordinates ...")
        # Main loop of creation
        for t in range(d[2]):
            for x in range(d[0]):
                for y in range(d[1]):

                    # Retreive voxel values
                    dx = voxels[x, y, t, 0]
                    dy = voxels[x, y, t, 1]
                    dt = voxels[x, y, t, 2] 
                    tb = values[x, y, t]

                    # STEP 1 - Convert transverse comoving distances to flat angular resolution
                    dz_pix = dz_val+dt/2 # Alter it by the CENTRAL pixel value

                    # Calculate transformed dimensions
                    dtheta = np.arctan(dx/dz_pix)
                    dphi = np.arctan(dy/dz_pix)

                    # STEPS 2 & 3 - Convert line-of-sight comoving distance to frequency
                    df = 0

                    if x == y == 0:
                        # Determine corresponding redshifts
                        z_bot = z_prev
                        z_top = cosmology.dz_to_z(dz+dt*u.Mpc)

                        # Convert to frequency
                        f_bot = cosmology.z_to_f(z_bot)
                        f_top = cosmology.z_to_f(z_top)

                        # Store altered frequency bandwidth
                        df = np.abs((f_bot-f_top).to_value(u.Hz))

                        # STEP 7.1 - Check if Regridding is needed
                        if df > max_freq_res_hz:
                            regrid_flag = regrid_flag or True

                    else:
                        df = voxels[0, 0, t, 2]

                    # STEP 4 - Convert brightness temperature to pixel flux
                    fxy = fq + df/2
                      
                    # Calculate pixelated luminosity
                    fv = SimulationReformatter.brightness_temperature_to_flux(tb=tb, fxy=fxy, dtheta=dtheta, dphi=dphi)
                
                    # Save values
                    voxels[x, y, t, 0] = dtheta
                    voxels[x, y, t, 1] = dphi
                    voxels[x, y, t, 2] = df
                    values[x, y, t] = fv
            
            # Increment distance and frequency parameters
            dz = dz + dt*u.Mpc # Increment the value of dz by voxel dimension
            z_prev = z_top # Top z in current box = bottom z in next box
            dz_val = dz_val + dt # Increment the value of dz by voxel dimension
            fq = fq + df # Increment cumulative frequency

            print("\rTransforming t-dim slice #", t, "of", d[2], end="")

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

                print("\rRegridding spaxel # (", x, ",", y, ")", "of", "(", d[0], ",", d[1], ")", end="")

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

                # Perform Regrid
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
    def calculate_cumulative_voxels(voxels, f_ref = 200 * u.MHz, phase_ref_point = omath.ZENITH_530):
        """
        Calculate the cumulative voxel sum and centre with a refrence point and frequency.

        :param voxels: An array describing a series of voxel dimensions corresponding to each sky model datacube voxel element in units of (rad, rad, Hz).
        :param f_ref: Refrence frequency, the ending frequency of the model.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.

        :return (rasum, decsum, freqsum): An array determining the specific central value of each voxel in its corresponding values array in units of (deg, deg, Hz).
        """

        # RA, Dec centering function and freq shift function
        def centering(x):
            """ Add half the total and half the spaxel widths to centre the main point. """
            return x - np.max(x, axis=(0, 1), keepdims=True)/2 - np.min(x, axis=(0, 1), keepdims=True)/2

        def shifting(x):
            """ Add only half the cell bandwidth, do not center. """
            return x - np.min(x, axis=2, keepdims=True)/2

        # Centre the cumulative RA and Dec sums so that the zero value is in the centre
        rasum = centering(np.cumsum(voxels[:,:,:,0], axis=0))
        decsum = centering(np.cumsum(voxels[:,:,:,1], axis=1))

        # Centre the frequency
        freqsum = f_ref.to_value(u.Hz) - shifting(np.cumsum(voxels[:,:,:,2], axis=2))

        # Calculate phase centre offsets
        source_pos = phase_ref_point.spherical_offsets_by(rasum * u.rad, decsum * u.rad)
        ras = source_pos.ra.to_value(u.deg)
        dcs = source_pos.dec.to_value(u.deg)

        return (ras, dcs, freqsum)

        
    @staticmethod
    # TODO: Gaussian source?
    def save_datacube_to_osm(values, voxels = None, cumulative_voxels = None, sigma_f = None, f_ref = 200 * u.MHz, phase_ref_point = omath.ZENITH_530, osm_output="reformat/osm_output/osm_output.osm"):
        """
        Saves a given datacube of flux values and voxel dimensions (RA, Dec, Freq.) to a master OSM file.

        :param values: The sky model datacube.
        :param voxels: An array describing a series of voxel dimensions corresponding to each sky model datacube voxel element (rad, rad, Hz).
        :param cumulative_voxels: A jagged array consisting of the cumulative summation of components from the voxel array (deg, deg, Hz).
        :param sigma_f: An array of same dimensions as values but containing information about the linewidth of the frequency emission profile.
        :param f_ref: Refrence frequency, the ending frequency of the model.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always reformat frequency bins, if false, reformat only when max frequency resolution is met.
        :param osm_output: The relative path to save the osm file to.
        """

        print("Configuring datacube for OSKAR file format ...")

        # Configure d variable
        d = np.shape(values)

        # Configure OSM path
        osm_output = ofc.expand_path(osm_output)

        # Cumulative sums are more important than voxel bins now
        (ras, dcs, freqsum) = (None, None, None) # Keep Pylint Happy
        if cumulative_voxels is None and voxels is None:
            raise ValueError("Either an array of voxels or cumulative voxes must be provided!")
        elif voxels is None:
            (ras, dcs, freqsum) = cumulative_voxels
        elif cumulative_voxels is None:
            (ras, dcs, freqsum) = SimulationReformatter.calculate_cumulative_voxels(voxels=voxels, f_ref=f_ref, phase_ref_point=phase_ref_point)
            
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
                        rascn = np.char.zfill(np.format_float_positional(ras[x, y, t], 6, False), 10)
                        decln = np.char.zfill(np.format_float_positional(np.abs(dcs[x, y, t]), 6, False), 9)
                        value = np.format_float_scientific(values[x, y, t], 4, False)
                        freq0 = np.format_float_positional(freqsum[x, y, t] / 1e6, 3, False)
                        linew = np.format_float_scientific(sigma_f[x, y, t], 4, False)

                        # Add +/- value to Declinations
                        if dcs[x, y, t] >= 0: decln = "+" + str(decln)
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

                    print("\rSaving data for spaxel # (", x, ",", y, ")", "of", "(", d[0], ",", d[1], ")", end="")

        print("\nProcess complete, data saved to "+osm_output)

    # TODO: Automatically find ideal UTC time of observation
    # pylint: disable=unused-argument
    @staticmethod
    def calculate_observation_time_from_date(phase_ref_point = omath.ZENITH_530, ref_time = omath.REF_TIME, ref_location = omath.SKA_REF_LOC, observation_length = omath.OBS_LEN_4HR):
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
    def generate_dynamic_settings(values, voxels = None, cumulative_voxels = None, phase_ref_point = omath.ZENITH_530, f_ref = 200 * u.MHz, ref_time = omath.REF_TIME, ref_location = omath.SKA_REF_LOC, observation_length = omath.OBS_LEN_4HR, save_dynamic_settings = ""):
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

        :return dynamic_settings: The dynamically defined settings dictionary.
        """
        # Calculated settings: [observation] start_frequency_hz, num_channels, frequency_inc_hz, phase_centre_ra_deg, phase_centre_dec_deg, length, start_time_utc
        # Calculated settings: [image] fov_deg, size

        # SETUP
        # Cumulative sums are more important than voxel bins now
        (ras, dcs, freqsum) = (None, None, None) # Keep Pylint Happy
        if cumulative_voxels is None and voxels is None:
            raise ValueError("Either an array of voxels or cumulative voxes must be provided!")
        elif voxels is None:
            (ras, dcs, freqsum) = cumulative_voxels
        elif cumulative_voxels is None:
            (ras, dcs, freqsum) = SimulationReformatter.calculate_cumulative_voxels(voxels=voxels, f_ref=f_ref, phase_ref_point=phase_ref_point)

        # Configure d variable
        d = np.shape(values)

        # Create deep copy of union/logical or settings set
        dynamic_settings = dict(ohelp.DEFAULT_GENERAL_SETTINGS)

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
        ref_time, _ = SimulationReformatter.calculate_observation_time_from_date(phase_ref_point=phase_ref_point, ref_time=ref_time, ref_location=ref_location, observation_length=observation_length)

        # Calculate RA dimension - use circular mean for angular averages
        rac = omath.diff(omath.circmean(ras[-1,:,:]), omath.circmean(ras[0,:,:]))

        # Calculate Dec dimension - use circular mean for angular averages
        decc = omath.diff(omath.circmean(dcs[-1,:,:], (-180, 180)), omath.circmean(dcs[0,:,:], (-180, 180)))

        # Set the field of view
        dynamic_settings['image']['fov_deg'] = max(rac, decc)

        # Set the observation time and length
        dynamic_settings['observation']['start_time_utc'] = (ref_time - observation_length / 2).utc.value
        dynamic_settings['observation']['length'] = str(observation_length.to_value(format='datetime'))

        # Save the dynamic settings
        if save_dynamic_settings != "":

            ofc.save_settings_from_dictionary(save_dynamic_settings, dynamic_settings)

            print("Saved dynamic and default settings to ini file: "+save_dynamic_settings)
        
        return dynamic_settings

    @staticmethod
    def generate_osm_from_simulation(values, voxels = None, z_ref = 7, phase_ref_point = omath.ZENITH_530, require_regrid = True, max_freq_res = 100 * u.MHz, v = (1.5, 1.5, 1.5), osm_output="reformat/osm_output/osm_output.osm", cosmology=eorcosmo(), save_dynamic_settings = "", ref_time = omath.REF_TIME, ref_location = omath.SKA_REF_LOC, observation_length = omath.OBS_LEN_4HR):
        """
        Generate a set of .osm files for an OSKAR sky model based on a Mpc**3 simulation output.

        :param values: The simulation datacube.
        :param z_ref: Refrence redshift, the ending redshift of the simulation.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always reformat frequency bins, if false, reformat only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution.
        :param v: If all voxels are the same, provides the initial voxel dimensions in h^-1 Mpc in dimensions (x, y, t), and auto-generates the voxel configuration array.
        :param osm_output: The relative path to save the osm file to.
        :param cosmology: The specific cosmology parameters in the form of a custom eorcosmo object.
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
        values, voxels, sigma_f, f_ref, regrid_flag = SimulationReformatter.transform_datacube_units(values=values, voxels=voxels, z_ref=z_ref, require_regrid=require_regrid, max_freq_res=max_freq_res, cosmology=cosmology)

        # STEP 7 - Regrid frequency-dimension data if needed
        if regrid_flag:
            values, voxels, sigma_f = SimulationReformatter.regrid_datacube(values=values, voxels=voxels, sigma_f=sigma_f, max_freq_res=max_freq_res)
        else:
            print("No Regrid required!")

        # STEP 8 - Write data to OSM file
        if osm_output != "":
            SimulationReformatter.save_datacube_to_osm(values=values, voxels=voxels, sigma_f=sigma_f, f_ref=f_ref, phase_ref_point=phase_ref_point, osm_output=osm_output)

        # Output dynamic settings file
        dynamic_settings = SimulationReformatter.generate_dynamic_settings(
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
    def generate_osm_from_h5(file, phase_ref_point = omath.ZENITH_530, require_regrid = True, max_freq_res = 100e6, osm_output="reformat/osm_output/osm_output.osm", coeval=True, ref_time = omath.REF_TIME, ref_location = omath.SKA_REF_LOC, observation_length = omath.OBS_LEN_4HR, save_dynamic_settings = ""):
        """
        Combines both the convert_h5_to_csv and generate_osm_from_simulation functions.

        :param file: Location of the h5 file.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always reformat frequency bins, if false, reformat only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution in Hz.
        :param osm_output: The directory to output the osm file.
        :param ref_time: An astropy.time.Time object stating the desired mid-observation time.
        :param ref_location: An astropy.coordinates.EarthLocation object stating the location of the telescope on Earth.
        :param observation_length: An astropy.time.TimeDelta object that gives the length of the observation.
        :param save_dynamic_settings: If non-empty, save the dynamic settings to an .ini file given by the path entered.

        :return dynamic_settings: The dynamically defined settings dictionary.
        """

        values, z_ref, vox, cosmology = SimulationReformatter.convert_h5_to_csv(file, coeval=coeval)

        if osm_output == "": osm_output = file.split('/')[-1][:-3] + "_osm.osm"

        return SimulationReformatter.generate_osm_from_simulation(values,
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
    def convert_osm_file_to_arrays(osm_file, generate_dynamic_settings = True, phase_ref_point_override = None, ref_time = omath.REF_TIME, ref_location = omath.SKA_REF_LOC, observation_length = omath.OBS_LEN_4HR, save_dynamic_settings = "", d = None):
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

        :return output_data: The values, voxels, cumulative voxels, sigma_f, f_ref, and phase_ref_point contained in a dictionary.
        
        If generate_dynamic_settings is set to True, additionally return an updated dictionary of settings to provide to OSKAR.
        """

        df = pd.read_csv(osm_file, delimiter=" ", skiprows=3, index_col=False, names=["RA", "Dec", "Stokes I", "Q", "U", "V", "Freq0"])

        #values, voxels = None, cumulative_voxels = None, phase_ref_point = omath.ZENITH_530, f_ref = 200 * u.MHz
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
            output_data["phase_ref_point"] = SkyCoord(ra=omath.denorm(omath.circmean(temp_values[:,:,:,0]))*u.deg, dec=omath.norm(omath.circmean(temp_values[:,:,:,1], (-180, 180)))*u.deg, frame='icrs') 
        else:
            output_data["phase_ref_point"] = phase_ref_point_override

        # Calculate the frequency bin step
        step = np.mean(temp_values[:,:,0, 2] - temp_values[:,:,1, 2])/2

        # Calculate the refrence frequency
        output_data["f_ref"] = np.mean(temp_values[:,:,0, 2]) + step

        # Begin mutating temp_values

        # Shift values based on refrence points
        temp_values[:,:,:,0] = omath.diff_sgn(temp_values[:,:,:,0], output_data["phase_ref_point"].ra.deg)
        temp_values[:,:,:,1] = omath.diff_sgn(temp_values[:,:,:,1], output_data["phase_ref_point"].dec.deg)
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
            dynamic_settings = SimulationReformatter.generate_dynamic_settings(
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
