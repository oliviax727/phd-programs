#!.venv/bin/python
"""
Run CHIPS.
"""

import sys

# sys.path.append("/home/olivia/Desktop/Programs/phd-programs/")
sys.path.append("/scratch/mwaeor/ohrw/phd-programs/")
from oskareor.power_spectra_analysis import ChipsManager as cmgr

SETONIX_DATA_DIR_SOFTWARE = "/software/projects/mwaeor/ohrw/"
SETONIX_DATA_DIR_SCRATCH = "/scratch/mwaeor/ohrw/"

# Locations

uvfits_loc = SETONIX_DATA_DIR_SCRATCH + "oskareor.data/oskar.data/fiducial/fiducial_uvw_plane.uvfits"
temp_parent_dir = SETONIX_DATA_DIR_SCRATCH
outdir = SETONIX_DATA_DIR_SCRATCH + "oskareor.data/oskar.data/fiducial/ps2ds/"

# Arguments
args = sys.argv

print("Arguments detected: " + str(args))

match int(args[1]):
    # Full CHIPS program
    case 0:
        print("Standard option selected")
        cmgr.run_chips_on_model(uvfits_loc, temp_parent_dir, outdir)

    # Chunk a file
    case 1:
        sim = args[2]
        print(f"Simulation {sim} selected")
        uvfits_loc = SETONIX_DATA_DIR_SCRATCH + f"oskareor.data/oskar.data/{sim}/{sim}_uvw_plane.uvfits"

        print("Chunking option selected")
        cmgr.extract_uvfits_chunks(uvfits_loc, temp_parent_dir)

    # Extract a file
    case 2:
        sim = args[2]
        print(f"Simulation {sim} selected")
        outdir = SETONIX_DATA_DIR_SCRATCH + f"oskareor.data/oskar.data/{sim}/ps2ds/"

        print("Extraction option selected")

        [f_chn, f_ref, index] = args[3:6]

        if len(args) > 6:
            nchan = args[6]
            kperp = args[7]
        else:
            nchan = kperp = 64

        cmgr.extract_and_save_2dps(
            outdir + f"{sim}_ps2d_data_chunk_{int(index)}.npz",
            temp_parent_dir=temp_parent_dir,
            f_ref=float(f_ref),
            f_chn=float(f_chn),
            n_channels=int(nchan),
            k_perp_bins=int(kperp),
        )

    # Chunk a template
    case 3:
        template = args[2]
        print(f"Template {template} selected")
        uvfits_loc = SETONIX_DATA_DIR_SCRATCH + f"oskareor.data/uvfits_templates/{template}_uvw_plane.uvfits"

        print("Chunking option selected")
        cmgr.extract_uvfits_chunks(uvfits_loc, temp_parent_dir, file_is_chunk=True)

    # Extract a template
    case 4:
        template = args[2]
        print(f"Template {template} selected")
        outdir = SETONIX_DATA_DIR_SCRATCH + f"oskareor.data/2dps_npz_templates/{template}_ps2ds/"

        print("Extraction option selected")

        [f_chn, f_ref, index] = args[3:6]

        if len(args) > 6:
            nchan = args[6]
            kperp = args[7]
        else:
            nchan = kperp = 64

        cmgr.extract_and_save_2dps(
            outdir + f"{template}_ps2d_data_chunk_{int(index)}.npz",
            temp_parent_dir=temp_parent_dir,
            f_ref=float(f_ref),
            f_chn=float(f_chn),
            n_channels=int(nchan),
            k_perp_bins=int(kperp),
        )

    # No option
    case -1:
        exit
