"""
The oskar_exec module contains all funtions that help to automate the execution of OSKAR - converting sky models to measurement sets and images for the analysis of simulated SKA observations.
"""

# System imports
import os
import subprocess

# Local imports
from oskar_helpers import OSKARHelper
from reformatter import Reformat

# TODO: Turn into pip project (later)
# TODO: Pydoc types and return values

# FIXME: Test refactored code

class BTAnalysisPipeline(object):
    """
    A broader class that combines all components of the individual components of the simulated IGM to simulated observation pipeline together.
    """

    @staticmethod
    def configure_oskar_settings(dynamic_settings = OSKARHelper.DEFAULT_GENERAL_SETTINGS, interferometer_settings_override = "", imager_settings_override = "", save_ini=""):
        """
        Configure the settings files for the OSKAR interferometer and imager programs.

        :param osm_file: File path to the existing OSM file
        :param dynamic_settings: Dynamically generated settings from the reformatting code.
        :param interferometer_settings_override: The file location of the OSKAR imager settings file. Leave blank if no override.
        :param imager_settings_override: The file location of the OSKAR interferometer settings template file. Leave blank if no override.
        :param use_imager: Whether or not to generate a dirty image with oskar_imager.
        :param save_ini: File path to save the compiled ini file. If blank, pass the settings only as a return.

        :return: The updated settings dictionary.
        """

        # Function to mutate the existing dynamic settings dictionary
        def mutate_settings(override_file, settings_dict):
            override_data = OSKARHelper.read_settings_to_dictionary(override_file)

            return settings_dict | override_data | OSKARHelper.PRIMARY_GENERAL_SETTINGS

        # Setup the interferometer ini file
        if interferometer_settings_override != "":
            dynamic_settings = mutate_settings(interferometer_settings_override, dynamic_settings)
            
        if imager_settings_override != "":
            dynamic_settings = mutate_settings(imager_settings_override, dynamic_settings)

        if save_ini != "":
            OSKARHelper.save_settings_from_dictionary(save_ini, dynamic_settings)
        
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
            settings[1] = OSKARHelper.read_settings_to_dictionary(settings[0])

        # Make sure that the settings exists before splitting
        if not settings[1] is None:
            interf_settings_dict = {k: settings[k] for k in OSKARHelper.LEGAL_INTERFEROMETER_HEADINGS}

            if use_imager:
                imager_settings_dict = {k: settings[k] for k in OSKARHelper.LEGAL_IMAGER_HEADINGS}

        # Make sure a file path has been provided to save the file to
        if settings[0] != "" and save_file:
            interf_settings_path = settings[0] + ".oskar_sim_interferometer.ini"

            OSKARHelper.save_settings_from_dictionary(interf_settings_path, interf_settings_dict)

            if use_imager:
                imager_settings_path = settings[0] + ".oskar_imager.ini"
                
                OSKARHelper.save_settings_from_dictionary(imager_settings_path, imager_settings_dict)

        return (interf_settings_path, interf_settings_dict), (imager_settings_path, imager_settings_dict)

    @staticmethod
    def run_oskar_on_osms(osm_file, interferometer_settings = ("", OSKARHelper.DEFAULT_INTERFEROMETER_SETTINGS), imager_settings = ("", OSKARHelper.DEFAULT_IMAGER_SETTINGS), fits_output="./fits_output.fits", oskar_exec=None, oskar_mode="python", use_imager=True):
        """
        Run oskar on each of the OSM sky models found in a fits directory, should already be formatted according to the output of the Reformat object.

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
                subprocess.run(["singularity","exec","--nv","--bind",cwd,"--cleanenv","--home",cwd,oskar_exec,"oskar_sim_interferometer",interferometer_settings[0]], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
            elif oskar_mode == "binary":
                subprocess.run([oskar_exec+"/oskar_sim_interferometer",interferometer_settings[0]], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
        except subprocess.CalledProcessError as e:
            print(f"Command '{e.cmd}' returned non-zero exit status {e.returncode}.")
            print(f"Error output: {e.stderr.decode()}")

        # Run OSKAR's imager simulation
        if use_imager:
            try:
                print("Running imager on "+osm_file)
                if oskar_mode == "singularity":
                    subprocess.run(["singularity","exec","--nv","--bind",cwd,"--cleanenv","--home",cwd,oskar_exec,"oskar_imager",imager_settings[0]], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
                elif oskar_mode == "binary":
                    subprocess.run([oskar_exec+"/oskar_imager",imager_settings[0]], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
            except subprocess.CalledProcessError as e:
                print(f"Command '{e.cmd}' returned non-zero exit status {e.returncode}.")
                print(f"Error output: {e.stderr.decode()}")

        subprocess.run(["find",".","-name","'*.log'","-type","f","-delete"], check=True, cwd=cwd)
        subprocess.run(["cp","BTA/oskar_output/sim_image_I.fits",fits_output], check=True)

    @staticmethod
    def setup_bta_dir(h5_file, interferometer_settings_override="./reformat/test_intif_inis/test_img_gen.ini", imager_settings_override="./reformat/test_intif_inis/test_intif_gen.ini", oskar_telescope_model="./oskar_run_stage/telescope_model_AAstar", template=False):
        """
        Sets up the operating directory from which all anaysis will be done.

        :param h5_file: The location of the h5 file.
        :param cd_in: Whether to cd into the directory once finished or not.
        :param interferometer_settings_override: The file location of the OSKAR imager settings file. Leave blank if no override.
        :param imager_settings_override: The file location of the OSKAR interferometer settings template file. Leave blank if no override.
        :param oskar_telescope_model: The telescope model for OSKAR to use.
        :param template: If true, handle and return no h5 data.
        :return: The new location of the h5, imager ini, and interterometer template ini files, as well as the telescope model location. Also returns the PWD if cd_in is True.
        """

        cwd = os.getcwd()

        # 1. Create Directory
        subprocess.run(["mkdir","-p","BTA"], check=True)

        # 2. Move h5 file and INIs to directory
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
    def h5_box_to_datacube(file, phase_ref_point = OSKARHelper.ZENITH_530, require_regrid = True, max_freq_res = 100e6, interferometer_settings_override = "./reformat/test_intif_inis/test_img_gen.ini", imager_settings_override = "./reformat/test_intif_inis/test_intif_gen.ini", outdir = ".", clean = True, oskar_exec = OSKARHelper.OSKAR_SIF, oskar_mode="singularity", oskar_telescope_model = OSKARHelper.TELESCOPE, template_preset = "", coeval = True, load_osm=False, ref_time = OSKARHelper.REF_TIME, ref_location = OSKARHelper.SKA_REF_LOC, observation_length = OSKARHelper.OBS_LEN_4HR, use_imager = True):
        """
        Full pipeline function for transforming a h5 simulation box output into a FITS datacube.

        :param file: Location of the h5 file.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always reformat frequency bins, if false, reformat only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution in Hz.
        :param interferometer_settings_override: The file location of the OSKAR imager settings file.
        :param imager_settings_override: The file location of the OSKAR interferometer settings template file.
        :param outdir: The output location of the final FITS file.
        :param oskar_exec: The SIF file or location of compiled OSKAR binaries.
        :param oskar_mode: How shall OSKAR be run? Options include: python, binary, command, singularity.
        :param oskar_telescope_model: The telescope model for OSKAR to use.
        :param template_preset: Use a mock values array instead of a h5 file. Ignores any provided h5 file.
        :param coeval: If the h5 box is coeval or lightcone based.
        :param load_osm: If true load treat the file variable as if it were an OSM file.
        :param ref_time: An astropy.time.Time object stating the desired mid-observation time.
        :param ref_location: An astropy.coordinates.EarthLocation object stating the location of the telescope on Earth.
        :param observation_length: An astropy.time.TimeDelta object that gives the length of the observation.
        :param use_imager: Whether or not to generate a dirty image with oskar_imager.
        """

        # Expand paths
        if oskar_mode == "binary" or oskar_mode == "singularity":
            oskar_exec = OSKARHelper.expand_path(oskar_exec)
        
        file = OSKARHelper.expand_path(file)
        imager_settings_override = OSKARHelper.expand_path(imager_settings_override)
        interferometer_settings_override = OSKARHelper.expand_path(interferometer_settings_override)
        oskar_telescope_model = OSKARHelper.expand_path(oskar_telescope_model)
        outdir = OSKARHelper.expand_path(outdir)

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
        osm_output = OSKARHelper.expand_path("BTA/" + h5_id + "_sky_model.osm")
        ini_output = OSKARHelper.expand_path("BTA/" + h5_id + "_general_settings.ini")
        fits_output = OSKARHelper.expand_path("BTA/" + h5_id + "_datacube.fits")

        # Set the default dynamic settings array
        dynamic_settings = OSKARHelper.DEFAULT_GENERAL_SETTINGS
        
        # Run the OSM and Settings generators
        if not os.path.isfile(osm_output):
            if not load_osm:
                if template_flag:
                    # IF we want to generate a fresh osm file AND its from a template
                    print("Generating OSM files from template ...")
                    template_values = Reformat.mock_values(template_preset)
                    dynamic_settings = Reformat.generate_osm_from_simulation(
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
                    print("Generating OSM files from h5 ...")
                    dynamic_settings = Reformat.generate_osm_from_h5(
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
                    subprocess.run(["cp", OSKARHelper.expand_path("~/.oskar/osm_templates/"+template_preset+"_sky_model.osm"), osm_output], check=True)
                    subprocess.run(["cp", OSKARHelper.expand_path("~/.oskar/ini_templates/"+template_preset+"_general_settings.ini"), ini_output], check=True)

                    dynamic_settings = OSKARHelper.read_settings_to_dictionary(ini_output)
                else:
                    # IF we want to skip generating the osm file AND a use a specified already-complete osm file
                    subprocess.run(["cp", file, osm_output], check=True)

                    _, dynamic_settings = Reformat.convert_osm_file_to_arrays(
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

class LoadDefaults:
    """
    A module for refreshing default sky models and measurement sets.
    """

    @staticmethod
    def load_default_sky_models():
        """
        (Re)load all default sky models and update corresponding ini and osm files.
        """

        for template_presett in Reformat.TEMPLATE_PRESETS:
            if "coeval" in template_presett:
                template_value = Reformat.mock_values(template_presett, d=(400, 400, 400))
            else:
                template_value = Reformat.mock_values(template_presett, scale=20)

            Reformat.generate_osm_from_simulation(
                template_value,
                osm_output=OSKARHelper.expand_path("~/.oskar/osm_templates/"+template_presett+"_sky_model.osm"),
                save_dynamic_settings=OSKARHelper.expand_path("~/.oskar/ini_templates/"+template_presett+"_general_settings.ini")
                )
    
    @staticmethod
    def load_default_oskar():
        """
        (Re)load all default sky models and update corresponding fits image and simulation .ms files.
        """

        for template_presett in Reformat.TEMPLATE_PRESETS:
            if "coeval" in template_presett:
                template_value = Reformat.mock_values(template_presett, d=(400, 400, 400))
            else:
                template_value = Reformat.mock_values(template_presett, scale=20)

            Reformat.generate_osm_from_simulation(
                template_value,
                osm_output=OSKARHelper.expand_path("~/.oskar/osm_templates/"+template_presett+"_sky_model.osm"),
                save_dynamic_settings=OSKARHelper.expand_path("~/.oskar/ini_templates/"+template_presett+"_general_settings.ini")
                )
