"""
The oskar_exec module contains all funtions that help to automate the execution of OSKAR - converting sky models to measurement sets and images for the analysis of simulated SKA observations.
"""

# System imports
import os
import subprocess

# Mathematics and calculations
import astropy.units as u
import numpy as np

# Astropy Extras
from astropy.units import Quantity

# Local imports
from oskareor.skalow_calc import OSKARFileConfig as ofc, SKAMath as omath
from oskareor.oskar_helpers import OSKARHelper as ohelp
from oskareor.reformatter import SimulationReformatter as simref

# TODO: Turn into pip project (later)
# TODO: Turn large parameter sets into more compartmentalised dictionaries (as suggested by PyLint)
# TODO: Build parameter error guards/handling for each function

# FIXME: Test refactored code

class BTAnalysisPipeline():
    """
    A broader class that combines all components of the individual components of the simulated IGM to simulated observation pipeline together.
    """

    @staticmethod
    def configure_oskar_settings(dynamic_settings = ohelp.DEFAULT_GENERAL_SETTINGS, interferometer_settings_override = "", imager_settings_override = "", save_ini=""):
        """
        Configure the settings files for the OSKAR interferometer and imager programs.

        :param osm_file: File path to the existing OSM file
        :param dynamic_settings: Dynamically generated settings from the reformatting code.
        :param interferometer_settings_override: The file location of the OSKAR imager settings file. Leave blank if no override.
        :param imager_settings_override: The file location of the OSKAR interferometer settings template file. Leave blank if no override.
        :param use_imager: Whether or not to generate a dirty image with oskar_imager.
        :param save_ini: File path to save the compiled ini file. If blank, pass the settings only as a return.

        :return dynamic_settings: The updated settings dictionary.
        """

        # Function to mutate the existing dynamic settings dictionary
        def mutate_settings(override_file, settings_dict):
            override_data = ofc.read_settings_to_dictionary(override_file)

            return settings_dict | override_data | ohelp.PRIMARY_GENERAL_SETTINGS

        # Setup the interferometer ini file
        if interferometer_settings_override != "":
            dynamic_settings = mutate_settings(interferometer_settings_override, dynamic_settings)
            
        if imager_settings_override != "":
            dynamic_settings = mutate_settings(imager_settings_override, dynamic_settings)

        if save_ini != "":
            ofc.save_settings_from_dictionary(save_ini, dynamic_settings)
        
        return dynamic_settings
    
    @staticmethod
    def split_general_settings(settings = ("", None), use_imager = True, save_file=True):
        """
        Split a settings dictionary into its components for the OSKAR interferometer and OSKAR imager.

        :param settings: A tuple containing a filepath and dictionary, the code will always prioritise using an evaluation the dictionary if available.
        :param use_imager: Whether to create a seperate imager file path and settings dictionary pair.
        :param save_file: If true, save to a file with a path and base name given by settings.

        :return (interferometer_settings, imager_settings): The two tuples containing the split settings (settings file location, settings dictionary)
        """

        # Retrieve sub-elements
        (gen_file, gen_dict) = settings

        # Declare default return values
        interf_settings_dict = None
        imager_settings_dict = None

        interf_settings_path = ""
        imager_settings_path = ""

        # If file is provided but not a dictionary, read the file
        if gen_dict is None and gen_file != "":
            gen_dict = ofc.read_settings_to_dictionary(gen_file)

        # Make sure that the settings exists before splitting
        if not gen_dict is None:
            interf_settings_dict = {k: gen_dict[k] for k in ohelp.LEGAL_INTERFEROMETER_HEADINGS}

            if use_imager:
                imager_settings_dict = {k: gen_dict[k] for k in ohelp.LEGAL_IMAGER_HEADINGS}

        # Remove ini file name, keep containing directory path and base preset name
        gen_file = gen_file[:-21]

        # Make sure a file path has been provided to save the file to
        if gen_file != "" and save_file:
            interf_settings_path = gen_file + "_oskar_sim_interferometer.ini"

            ofc.save_settings_from_dictionary(interf_settings_path, interf_settings_dict)

            if use_imager:
                imager_settings_path = gen_file + "_oskar_imager.ini"
                
                ofc.save_settings_from_dictionary(imager_settings_path, imager_settings_dict)

        return (interf_settings_path, interf_settings_dict), (imager_settings_path, imager_settings_dict)

    @staticmethod
    def run_oskar_on_osm(osm_file, interferometer_settings = ("", ohelp.DEFAULT_INTERFEROMETER_SETTINGS), imager_settings = ("", ohelp.DEFAULT_IMAGER_SETTINGS), oskar_exec=None, oskar_mode="python", use_imager=True):
        """
        Run oskar on each of the OSM sky models found in a fits directory, should already be formatted according to the output of the SimulationReformatter object.

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

        cwd = os.getcwd()

        # Set up settings path/dict
        settings = [None, None]

        if oskar_mode != "python":
            # Expand path if not python
            settings[0] = ofc.expand_path(interferometer_settings[0])
            oskar_exec = ofc.expand_path(oskar_exec)

            if use_imager:
                settings[1] = ofc.expand_path(imager_settings[0])

        else:
            # Copy dictionary if python
            settings[0] = interferometer_settings[1]

            if use_imager:
                settings[1] = imager_settings[1]
            
        settings = tuple(settings)

        # Executes an shell command with all try-catches included
        def execute_oskar_shell_command(command):
            shell_command = " ".join(command)

            print(f"Attempting to run the command - $ {shell_command}")

            try:
                subprocess.run(shell_command, check=True, shell=True, stderr=subprocess.PIPE) #stdout=subprocess.PIPE
            except subprocess.CalledProcessError as e:
                print(f"Command '{e.cmd}' returned non-zero exit status {e.returncode}.")
                print(f"Error output: {e.stderr.decode()}")

        # Run OSKAR's interferometer simulation
        print("Running interferometer on "+osm_file)
        if oskar_mode == "singularity":
            execute_oskar_shell_command(["singularity","exec","--nv","--bind",cwd,"--cleanenv","--home",cwd,oskar_exec,"oskar_sim_interferometer",settings[0]])
        elif oskar_mode == "binary":
            execute_oskar_shell_command([oskar_exec+"/oskar_sim_interferometer",settings[0]])
        
        # Run OSKAR's imager simulation
        if use_imager:
            print("Running imager on "+osm_file)
            if oskar_mode == "singularity":
                execute_oskar_shell_command(["singularity","exec","--nv","--bind",cwd,"--cleanenv","--home",cwd,oskar_exec,"oskar_imager",settings[1]])
            elif oskar_mode == "binary":
                execute_oskar_shell_command([oskar_exec+"/oskar_imager",settings[1]])

    @staticmethod
    def setup_bta_dir(oskar_telescope_model, h5_file="", interferometer_settings_override="", imager_settings_override="", template=False):
        """
        Sets up the operating directory from which all anaysis will be done.

        :param oskar_telescope_model: The telescope model for OSKAR to use.
        :param h5_file: The location of the h5 file.
        :param interferometer_settings_override: The file location of the OSKAR imager settings file. Leave blank if no override.
        :param imager_settings_override: The file location of the OSKAR interferometer settings template file. Leave blank if no override.
        :param template: If true, handle and return no h5 data.

        :return (h5_file, interf_override_ini, imager_override_ini, telescope_model, cwd): The new location of the h5, imager ini, and interterometer template ini files, as well as the telescope model location.
        """

        cwd = os.getcwd()

        # Define default return value
        default_return = ["", "", "", "BTA/telescope_model.tm", cwd]

        # Create directory, delete contents of existing BTA directory
        subprocess.run(["rm","-rf","BTA"], check=True)
        subprocess.run(["mkdir","-p","BTA"], check=True)

        # Move h5 file and INIs to directory
        if not template and h5_file != "":
            subprocess.run(["cp",ofc.expand_path(h5_file),"BTA/analysis.h5"], check=True)
            default_return[0] = "BTA/analysis.h5"

        if interferometer_settings_override != "":
            subprocess.run(["cp",ofc.expand_path(interferometer_settings_override),"BTA/interferometer_override.ini"], check=True)
            default_return[1] = "BTA/interferometer_override.ini"

        if imager_settings_override != "":
            subprocess.run(["cp",ofc.expand_path(imager_settings_override),"BTA/imager_override.ini"], check=True)
            default_return[2] = "BTA/imager_override.ini"

        subprocess.run(["cp","-r",ofc.expand_path(oskar_telescope_model),"BTA/telescope_model.tm"], check=True)

        return default_return

    @staticmethod
    def clean_bta_dir(outpath, use_imager = True, clean = True, settings = None):
        """
        Finishes and cleans up the mess created by the BTA class.

        :param outpath: A tuple containing the output path of the: sim.ms file, vis file, and image FITS file.
        :param use_imager: Whether or not copy any existing imager FITS file.
        :param clean: Whether or not to remove the BTA directory. Only works if cd_out is true.
        :param settings: The settings used to run OSKAR, so that the function knows where to find the outputs.
        """

        # Keep pylint happy
        if settings is None:
            settings = ohelp.DEFAULT_GENERAL_SETTINGS

        # Get original BTA location string
        original_locations = np.array([settings["interferometer"]["ms_filename"], settings["interferometer"]["oskar_vis_filename"], settings["image"]["root_path"]+"_I.fits"])
        original_locations = tuple(map(ofc.expand_path, original_locations))

        # Convert outpath string
        outpath = tuple(map(ofc.expand_path, outpath))

        # Wrap in try-catch just in case oskar stops working
        try:
            # Move sim.ms file to directory
            if outpath[0] != "":
                subprocess.run(["mv", original_locations[0], outpath[0]], check=True)

            # Move vis file to directory
            if outpath[1] != "":
                subprocess.run(["mv", original_locations[1], outpath[1]], check=True)

            # Move image fits file to directory
            if outpath[2] != "" and use_imager:
                subprocess.run(["mv", original_locations[0], outpath[2]], check=True)

        except subprocess.CalledProcessError as e:
            print(f"Command '{e.cmd}' returned non-zero exit status {e.returncode}.")
            print(f"Error output: {e.output}")
        
        finally:
            if clean:
                subprocess.run(["rm","-rf","BTA"], check=True)

    @staticmethod
    def run_oskar_on_model(file="", phase_ref_point = omath.ZENITH_530, require_regrid = True, max_freq_res: Quantity = 100e6 * u.MHz, interferometer_settings_override = "", imager_settings_override = "", outpath = ("","",""), clean = True, oskar_exec = "", oskar_mode="singularity", oskar_telescope_model = "", template_preset = "", coeval = True, load_osm=False, ref_time = omath.REF_TIME, ref_location = omath.SKA_REF_LOC, observation_length = omath.OBS_LEN_4HR, use_imager = True, oskar_parent_dir = "~"):
        """
        Full pipeline function for transforming a h5 simulation box output into a FITS datacube.

        :param file: Location of the h5 file.
        :param phase_ref_point: An astropy.coordinates.SkyCoord object stating the central sky refrence point.
        :param require_regrid: If true then always reformat frequency bins, if false, reformat only when max frequency resolution is met.
        :param max_freq_res: Maximum allowable voxel frequency resolution in Hz.
        :param interferometer_settings_override: The file location of the OSKAR imager settings file.
        :param imager_settings_override: The file location of the OSKAR interferometer settings template file.
        :param outpath: A tuple containing the output path of the: sim.ms file, vis file, and image FITS file.
        :param oskar_exec: The SIF file or location of compiled OSKAR binaries.
        :param oskar_mode: How shall OSKAR be run? Options include: python, binary, command, singularity.
        :param oskar_telescope_model: The telescope model for OSKAR to use.
        :param template_preset: Use a mock values array instead of a h5 file. Ignores any provided h5 file.
        :param coeval: If the h5 box is coeval or lightcone based.
        :param load_osm: If true treat the file variable as if it were an OSM file.
        :param ref_time: An astropy.time.Time object stating the desired mid-observation time.
        :param ref_location: An astropy.coordinates.EarthLocation object stating the location of the telescope on Earth.
        :param observation_length: An astropy.time.TimeDelta object that gives the length of the observation.
        :param use_imager: Whether or not to generate a dirty image with oskar_imager.
        :param oskar_parent_dir: The directory containing the .oskar folder (default is the home folder).
        """

        # Set defaults
        if oskar_exec == "": oskar_exec = oskar_parent_dir + ohelp.OSKAR_SIF
        if oskar_telescope_model == "": oskar_telescope_model = oskar_parent_dir + ohelp.TELESCOPE

        template_flag = (template_preset != "")
        
        # Set templates
        if load_osm and file != "": template_preset = file.split('/')[-1][:-4]

        print("Setting up BTA directory ...")
        h5_file, interf_override_ini, imager_override_ini, _, _ = BTAnalysisPipeline.setup_bta_dir(
            oskar_telescope_model=oskar_telescope_model,
            h5_file=file,
            interferometer_settings_override=interferometer_settings_override,
            imager_settings_override=imager_settings_override,
            template=template_flag
            )

        # Create output file locations
        h5_id = file.split('/')[-1][:-3] if not template_flag else template_preset
        osm_output = ofc.expand_path("BTA/sky_model.osm")
        ini_output = ofc.expand_path("BTA/" + h5_id + "_general_settings.ini")

        # Set the default dynamic settings array
        dynamic_settings = ohelp.DEFAULT_GENERAL_SETTINGS
        
        # Run the OSM and Settings generators
        if not os.path.isfile(osm_output):
            if not load_osm:
                if template_flag:
                    # IF we want to generate a fresh osm file AND its from a template
                    print("Generating OSM files from template ...")
                    template_values = simref.mock_values(template_preset, oskar_parent_dir=oskar_parent_dir)
                    dynamic_settings = simref.generate_osm_from_simulation(
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
                    dynamic_settings = simref.generate_osm_from_h5(
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
                    subprocess.run(["cp", ofc.expand_path(ohelp.default_template_path(template_preset=template_preset, oskar_parent_dir=oskar_parent_dir, file_type="osm")), osm_output], check=True)
                    subprocess.run(["cp", ofc.expand_path(ohelp.default_template_path(template_preset=template_preset, oskar_parent_dir=oskar_parent_dir, file_type="ini")), ini_output], check=True)

                    dynamic_settings = ofc.read_settings_to_dictionary(ini_output)
                else:
                    # IF we want to skip generating the osm file AND a use a specified already-complete osm file
                    subprocess.run(["cp", ofc.expand_path(file), osm_output], check=True)

                    _, dynamic_settings = simref.convert_osm_file_to_arrays(
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
        BTAnalysisPipeline.run_oskar_on_osm(
            osm_file=osm_output,
            interferometer_settings=interferometer_settings,
            imager_settings=imager_settings,
            oskar_exec=oskar_exec,
            oskar_mode=oskar_mode,
            use_imager=use_imager
            )

        # If clean is true remove all data relating to execution
        print("Cleaning up ...")
        BTAnalysisPipeline.clean_bta_dir(outpath=outpath, use_imager=use_imager, clean=clean, settings=dynamic_settings)


class LoadDefaults:
    """
    A module for refreshing and reloading templates. By default all relevant files for all relevant templates will be reloaded.
    The list of available templates are found in the `TEMPLATES` set and the list of available file outputs are found in the `FILETYPES` set.

    NB: Note that each fits datacube will take a large quantity of space on disk >10 GB.
    """

    # Update Settings
    TEMPLATES = set(simref.TEMPLATE_PRESETS.keys())
    FILETYPES = { "osm", "ini", "ms", "vis", "fits" }

    @staticmethod
    def reload_template_sky_models(update_which_files = None, update_which_templates = None, oskar_parent_dir = "~"):
        """
        (Re)load all default sky models and update corresponding ini and osm files.

        :param update_which_files: A set containing all of the file types to be updated, see `FILETYPES` for what is available (and used as a default).
        :param update_which_templates: A set containing all of the templates to be updated, see `TEMPLATES` for what is available (and used as a default).
        :param oskar_parent_dir: The directory containing the .oskar folder (default is the home folder).
        """

        # Fallback to defaults
        if update_which_files is None:
            update_which_files = LoadDefaults.FILETYPES
        if update_which_templates is None:
            update_which_templates = LoadDefaults.TEMPLATES

        # Set update_which parameters to be lower case
        update_which_files = ofc.recase_iterable(update_which_files)
        update_which_templates = ofc.recase_iterable(update_which_templates)

        # Loop through all selected templates
        for template_preset_loop in update_which_templates:
            if "coeval" in template_preset_loop:
                template_value = simref.mock_values(template_preset_loop, oskar_parent_dir=oskar_parent_dir)
            else:
                template_value = simref.mock_values(template_preset_loop, scale=20, oskar_parent_dir=oskar_parent_dir)

            simref.generate_osm_from_simulation(
                template_value,
                osm_output=ohelp.default_template_path(template_preset=template_preset_loop, oskar_parent_dir=oskar_parent_dir, file_type="osm") if "osm" in update_which_files else "",
                save_dynamic_settings=ohelp.default_template_path(template_preset=template_preset_loop, oskar_parent_dir=oskar_parent_dir, file_type="ini") if "ini" in update_which_files else ""
                )
    
    @staticmethod
    def reload_template_oskar_sims(start_from_scratch = False, update_which_files = None, update_which_templates = None, oskar_parent_dir = "~"):
        """
        (Re)load all default sky models and update corresponding measurement sets, visibility tables, and fits datacube files.

        :param start_from_scratch: If true, run the all templates from scratch, and automatically generate the osms and sky models. This will not update the sky models in the template folder.
        :param update_which_files: A set containing all of the file types to be updated, see `FILETYPES` for what is available (and used as a default).
        :param update_which_templates: A set containing all of the templates to be updated, see `TEMPLATES` for what is available (and used as a default).
        :param oskar_parent_dir: The directory containing the .oskar folder (default is the home folder).
        :param default_settings_override: A modified version of the default settings, for experimentation purposes.
        """

        # Fallback to defaults
        if update_which_files is None:
            update_which_files = LoadDefaults.FILETYPES
        if update_which_templates is None:
            update_which_templates = LoadDefaults.TEMPLATES

        # Set update_which parameters to be lower case
        update_which_files = ofc.recase_iterable(update_which_files)
        update_which_templates = ofc.recase_iterable(update_which_templates)

        # Loop through all selected templates
        for template_preset_loop in update_which_templates:
            BTAnalysisPipeline.run_oskar_on_model(
                template_preset=template_preset_loop,
                outpath=(
                    ohelp.default_template_path(template_preset=template_preset_loop, oskar_parent_dir=oskar_parent_dir, file_type="ms") if "ms" in update_which_files else "",
                    ohelp.default_template_path(template_preset=template_preset_loop, oskar_parent_dir=oskar_parent_dir, file_type="vis") if "vis" in update_which_files else "",
                    ohelp.default_template_path(template_preset=template_preset_loop, oskar_parent_dir=oskar_parent_dir, file_type="fits") if "fits" in update_which_files else ""
                ),
                oskar_mode="binary",
                oskar_exec=oskar_parent_dir+ohelp.OSKAR_BIN,
                use_imager=(".fits" in update_which_files),
                load_osm=(not start_from_scratch),
                oskar_parent_dir=oskar_parent_dir
                )
                
    @staticmethod
    def reload_all(update_which_files = None, update_which_templates = None, oskar_parent_dir = "~"):
        """
        (Re)load all default sky models and update the corresponding sky models, settings files, measurement sets, visibility tables, and fits datacubes.

        :param use_imager: Whether or not to also generate a datacube of the template file (note that each fits datacube will take a large quantity of space on disk >10 GB).
        :param update_which_files: A set containing all of the file types to be updated, see `FILETYPES` for what is available (and used as a default).
        :param update_which_templates: A set containing all of the templates to be updated, see `TEMPLATES` for what is available (and used as a default).
        :param oskar_parent_dir: The directory containing the .oskar folder (default is the home folder).
        :param default_settings_override: A modified version of the default settings, for experimentation purposes.
        """

        # Fallback to defaults
        if update_which_files is None:
            update_which_files = LoadDefaults.FILETYPES
        if update_which_templates is None:
            update_which_templates = LoadDefaults.TEMPLATES

        # Set update_which parameters to be lower case
        update_which_files = ofc.recase_iterable(update_which_files)
        update_which_templates = ofc.recase_iterable(update_which_templates)

        # Call previous loader functions
        if "ini" in update_which_files or "osm" in update_which_files:
            LoadDefaults.reload_template_sky_models(update_which_files=update_which_files, update_which_templates=update_which_templates, oskar_parent_dir=oskar_parent_dir)
        if "ms" in update_which_files or "vis" in update_which_files or "fits" in update_which_files:
            LoadDefaults.reload_template_oskar_sims(update_which_files=update_which_files, update_which_templates=update_which_templates, oskar_parent_dir=oskar_parent_dir)
