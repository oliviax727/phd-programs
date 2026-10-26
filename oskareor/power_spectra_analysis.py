"""The power spectra analysis python module handles any and all analysis of power spectra outputted by CHIPS. Including plotting, condensing to a 1D spectrum, and statistical analysis."""

# System imports
import os
import subprocess

# Mathematics and calculations
from astropy import units as u, constants as c
import numpy as np
from pyuvdata import UVData

# Astropy Extras
from astropy.units import Quantity
from chips_wrappers.ps_methods.chips_data import ChipsDataProducts
from chips_wrappers.setup.make_args import FakeArgs

# Local imports
from oskareor.skalow_calc import FileManager as ofmg, SKAMath as omath, SKAString as ostr, EoRCosmology as Cosmo


class Chunker:
    """Converts or calculates data into chunks."""

    # Print channel centres, widths, and start/end points and append to an accumulator list

    @staticmethod
    def chunk_array(data, voxels, cosmo=None, n_chunks=32, show_chunks=False, z_ref=5.5):
        """Split a large 3-dimensional datacube into a series of overlapping, evenly distributed, and fixed-size chunks for analysis in a power spectra.

        :param data: The 3-dimensional array to divide into chunks.
        :param voxels: The 4-dimensional voxel array to divide into chunks.
        :param n_chunks: The number of chunks to split the data in.
        :param show_chunks: Whether to show the chunk boundaries, defaults to False
        :param z_ref: The starting redshift of the box.
        :param cosmo: The EoRCosmology object for the box. If none, don't perform redshift calculators.

        :return (chunks, voxels, z_bins): Two 4-dimensional arrays containing each 3-dimensional chunk, one containing temperature values in Kelvin and the other voxel dimensions in cMpc. The third argument describes the start, middle, and end z-values for each chunk.
        """

        print("Dividing array into chunks ...")

        b = 2 * int(data.shape[2] / n_chunks)  # Chunk size
        dz_inc = np.mean(voxels[:, :, :, 2]) * u.Mpc  # Average voxel depth

        if cosmo is not None:
            dz_ref = cosmo.z_to_dz(z_ref) + dz_inc
            chunk_z_bins = []

        chunk_values = []
        chunk_voxels = []

        if show_chunks:
            print("Chunk data:")
            print("+" + "—" * 86 + "+")
            print("| Chunk No.        Range ([start, end])    Centre      Length    Redshift Bins         |")
            print("+" + "—" * 86 + "+")

        for k in range(n_chunks - 1):
            c = int(b / 2 * (k + 1))
            i = int(c - b / 2)
            f = int(c + b / 2)

            if cosmo is not None:
                z = np.array(
                    [
                        cosmo.dz_to_z(dz_ref + i * dz_inc),
                        cosmo.dz_to_z(dz_ref + c * dz_inc),
                        cosmo.dz_to_z(dz_ref + f * dz_inc),
                    ]
                )
            else:
                z = np.array([-1.0, -1.0, -1.0])

            if show_chunks:
                print(
                    f"| {k}"
                    + " " * (4 - len(str(k)))
                    + " " * 19
                    + f"[{i: 5d}, {f: 5d}]"
                    + " " * 5
                    + f"{c: 5d}"
                    + " " * 9
                    + f"{b:3d}    "
                    + f"[{z[0]:-5.2f}, {z[1]:-5.2f}, {z[2]:-5.2f}] |"
                )

            chunk_values.append(data[:, :, i:f])
            chunk_voxels.append(voxels[:, :, i:f])
            chunk_z_bins.append(z)

        chunk_values = np.array(chunk_values)
        chunk_voxels = np.array(chunk_voxels)
        if cosmo is not None:
            chunk_z_bins = np.array(chunk_z_bins)

        if show_chunks:
            print("+" + "—" * 86 + "+")

            print(f"Final post-chunk shape: {chunk_values.shape}")

        if cosmo is not None:
            return chunk_values, chunk_voxels, chunk_z_bins
        else:
            return chunk_values, chunk_voxels

    @staticmethod
    def generate_chunks(n_channels, voxel_depth=2, cosmo=None, n_chunks=32, show_chunks=False, z_ref=5.5):
        """A simpler method that does not split an array but instead returns the chunk bins in channels.

        :param n_channels: The mumber of frequency channels in the data array.
        :param n_chunks: The number of chunks to split the data in.
        :param voxel_depth: The depth of the voxels in cMpc
        :param show_chunks: Whether to show the chunk boundaries, defaults to False
        :param z_ref: The starting redshift of the box.
        :param cosmo: The EoRCosmology object for the box. If none, don't perform redshift calculators.

        :return (chunks, voxels, z_bins): Two 4-dimensional arrays containing each 3-dimensional chunk, one containing temperature values in Kelvin and the other voxel dimensions in cMpc. The third argument describes the start, middle, and end z-values for each chunk.
        """

        print("Dividing array into chunks ...")

        b = 2 * int(n_channels / n_chunks)  # Chunk size
        dz_inc = np.mean(voxel_depth) * u.Mpc  # Average voxel depth

        if cosmo is not None:
            dz_ref = cosmo.z_to_dz(z_ref) + dz_inc
            chunk_z_bins = []

        chunk_bounds = []

        if show_chunks:
            print("Chunk data:")
            print("+" + "—" * 86 + "+")
            print("| Chunk No.        Range ([start, end])    Centre      Length    Redshift Bins         |")
            print("+" + "—" * 86 + "+")

        for k in range(n_chunks - 1):
            c = int(b / 2 * (k + 1))
            i = int(c - b / 2)
            f = int(c + b / 2)

            if cosmo is not None:
                z = np.array(
                    [
                        cosmo.dz_to_z(dz_ref + i * dz_inc),
                        cosmo.dz_to_z(dz_ref + c * dz_inc),
                        cosmo.dz_to_z(dz_ref + f * dz_inc),
                    ]
                )
                chunk_z_bins.append(z)

            if show_chunks:
                print(
                    f"| {k}"
                    + " " * (4 - len(str(k)))
                    + " " * 19
                    + f"[{i: 5d}, {f: 5d}]"
                    + " " * 5
                    + f"{c: 5d}"
                    + " " * 9
                    + f"{b:3d}    "
                    + f"[{z[0]:-5.2f}, {z[1]:-5.2f}, {z[2]:-5.2f}] |"
                )

            chunk_bounds.append(np.array([i, c, f]))

        chunk_bounds = np.array(chunk_bounds)

        if cosmo is not None:
            chunk_z_bins = np.array(chunk_z_bins)

        if show_chunks:
            print("+" + "—" * 86 + "+")

        if cosmo is not None:
            return chunk_bounds, chunk_z_bins
        else:
            return chunk_bounds


class ChipsManager:
    """Runs and evaluates the CHIPS program"""

    DEFAULT_TEMP_FILE = os.getcwd() + "/PSA/chunk_data.uvfits"

    @staticmethod
    def extract_uvfits_chunk(uvfits_file, chunk_boundaries, outfile="./PSA/chunk_data.uvfits"):
        """Retrieves the UVFits file and then isolates the data to a specific chunk.

        .. NB: Channel cropping indexes from zero and isolates the data in a similar way to numpy array indexing, where the start channel is included but the end channel is not e.g. for `chunk_boundaries = (i, c, f)`, the saved data contains channels i to f-1 (inclusive to both ends).

        :param uvfits_file: The UVFits file location to convert.
        :param chunk_boundaries: A tuple containing the start, central, and end channels to crop to.
        """

        # Read the UV Fits file
        uvd = UVData()
        uvd.read_uvfits(ofmg.expand_path(uvfits_file))

        # Extract chunk data
        i, _, f = chunk_boundaries

        # Create a copy of the object
        uvd.select(freq_chans=np.arange(i, f))

        # Write the object to the new file
        uvd.write_uvfits(outfile)

    @staticmethod
    def run_chips(uvfits_file, temp_parent_dir=".", n_channels=64, k_perp_bins=64, u_max=1000.0):
        """Runs the CHIPS program, assumes the CHIPS directory has been created and the files are split.

        :param uvfits_file: The UVFits file location to run the full CHIPS program on.
        """

        os.environ["OUTPUTDIR"] = ofmg.expand_path(temp_parent_dir)
        os.environ["INPUTDIR"] = ofmg.expand_path(temp_parent_dir)

        subprocess.run(
            [
                "/software/projects/mwaeor/ctrott/setonix/chips_2025_ska/gridvisska",
                ofmg.expand_path(uvfits_file),
                "side",
                "15.0",  # Size of the primary beam in fourier space
                str(c.c.si.value),
            ],
            check=True,
        )

        subprocess.run(
            [
                "/software/projects/mwaeor/ctrott/setonix/chips_2025_ska/prepare_ska",
                "side",
                str(n_channels),
                "0",  # Indexing start value
                "'yy'",
                "side",
                "1",  # Band number - legacy from MWA
            ],
            check=True,
        )

        subprocess.run(
            [
                "/software/projects/mwaeor/ctrott/setonix/chips_2025_ska/lssa_fg_ska",
                "side",
                str(n_channels),
                str(k_perp_bins),
                "'yy'",
                str(u_max),
                "chips.out",
                "0",  # Bias mode, 0 means evaluate whole band
                "1",  # Band number - legacy from MWA
            ],
            check=True,
        )

    @staticmethod
    def extract_and_save_2dps(
        outdir,
        cosmo=None,
        temp_parent_dir=".",
        n_channels=64,
        k_perp_bins=64,
        u_max=1000,
        f_ref=220e6,
        f_chn=140e3,
    ):
        """Extracts the CHIPS output of a chunk, combines the data, and then saves it to a .npz file for later useage. Uses the `CHIPS_WRAPPERS` `ChipsDataProducts` object for posterity."""

        if cosmo is None:
            cosmo = Cosmo()

        args = FakeArgs()
        args.basedir = ofmg.expand_path(temp_parent_dir)
        args.chips_tag = "chips.out"
        args.polarisation = "YY"
        args.N_chan = n_channels
        args.N_kperp = k_perp_bins
        args.umax = u_max
        args.omega_matter = cosmo.omega_m_0
        args.omega_baryon = cosmo.omega_b_0
        args.omega_lambda = 1 - cosmo.omega_m_0
        args.hubble = cosmo.little_h * 100

        # Defaults
        args.lowerfreq = f_ref
        args.chan_width = f_chn
        args.Neta = args.N_chan // 2

        cdp = ChipsDataProducts(args)

        cdp._read_in_data_and_convert("yy", "chips.out", False, ofmg.expand_path(outdir))

    @staticmethod
    def run_chips_on_model(
        uvfits_file,
        temp_parent_dir=".",
        outdir=".",
        z_range=[5.5, 13.5],
        n_channels=1024,
        n_chunks=32,
        k_perp_bins=64,
        u_max=1000,
    ):
        """Runs the entirety of the CHIPS model."""

        print("Creating the output directory ./PSA ...")

        # Create the directory and/or clear it
        subprocess.run(["rm", "-rf", temp_parent_dir + "/PSA/chunks"], check=True)
        subprocess.run(["mkdir", "-p", temp_parent_dir + "/PSA/chunks"], check=True)

        print("Generating chunks ...")

        # Create an array of chunks
        chunk_bins = Chunker.generate_chunks(n_channels=n_channels, n_chunks=n_chunks, show_chunks=False)

        for i, chunk_bounds in enumerate(chunk_bins):
            print(f"Extracting and evaluating chunk: {i} of {n_chunks}")

            f_tot = Cosmo.z_to_f(z_range[1]) - Cosmo.z_to_f(z_range[0])
            f_chn = f_tot / n_channels
            f_ref = Cosmo.z_to_f(z_range[0]) - chunk_bounds[0] * f_chn

            # Create the temporary UVFits file location for CHIPS to analyse
            uvfits_loc = temp_parent_dir + f"/PSA/chunks/chunk{i}.uvfits"

            ChipsManager.extract_uvfits_chunk(
                uvfits_file=uvfits_file, chunk_boundaries=chunk_bounds, outfile=uvfits_loc
            )

            ChipsManager.run_chips(
                uvfits_loc,
                temp_parent_dir=temp_parent_dir + "/PSA/",
                n_channels=chunk_bounds[1] - chunk_bounds[0],
                k_perp_bins=k_perp_bins,
                u_max=u_max,
            )

            ChipsManager.extract_and_save_2dps(
                outdir + f"fiducial_ps2d_data_chunk_{i}.npz",
                temp_parent_dir=temp_parent_dir + "/PSA/",
                n_channels=chunk_bounds[1] - chunk_bounds[0],
                k_perp_bins=k_perp_bins,
                u_max=u_max,
                f_ref=f_ref,
                f_chn=f_chn,
            )

    def extract_uvfits_chunks(
        uvfits_file,
        temp_parent_dir=".",
        z_range=[5.5, 13.5],
        n_channels=1024,
        n_chunks=32,
    ):
        """Creates a series of UVFits chunk files to later evaluate."""

        print("Creating the output directory ./PSA ...")

        # Create the directory and/or clear it
        subprocess.run(["rm", "-rf", temp_parent_dir + "/PSA/chips_output"], check=True)
        subprocess.run(["mkdir", "-p", temp_parent_dir + "/PSA/chips_output"], check=True)

        print("Generating chunks ...")

        # Create an array of chunks
        chunk_bins = Chunker.generate_chunks(n_channels=n_channels, n_chunks=n_chunks, show_chunks=False)

        for i, chunk_bounds in enumerate(chunk_bins):
            print(f"Extracting and evaluating chunk: {i} of {n_chunks}")

            f_tot = Cosmo.z_to_f(z_range[1]) - Cosmo.z_to_f(z_range[0])
            f_chn = f_tot / n_channels
            f_ref = Cosmo.z_to_f(z_range[0]) - chunk_bounds[0] * f_chn

            # Create the temporary UVFits file location for CHIPS to analyse
            uvfits_loc = temp_parent_dir + f"/PSA/chunks/chunk_{f_chn}_{f_ref}_{i}_temp.uvfits"

            ChipsManager.extract_uvfits_chunk(
                uvfits_file=uvfits_file, chunk_boundaries=chunk_bounds, outfile=uvfits_loc
            )

    @staticmethod
    def extract_and_load_2dps(npz_file):
        """Reads chunk data saved to an `.npz` file and returns the array as if it was ran through `PSTools.cylindrical_power_spectrum`."""

        npz = npz = np.load(npz_file)

        return npz["twoD_power"] * npz["weights"], npz["k_perp"], npz["k_parr"]


class PSTools:
    """The PSTools class contains any and all methods to analyse the power spectra of data either pre- or post- observation."""

    @staticmethod
    def cubic_power_spectra(data, voxels, cosmo):
        """Calculate the cubic/3-dimensional power spectra for a given volume of data. Assumes that the voxel data is uniform.

        :param data: A 3-dimensional array of brightness temperature values in Kelvin.
        :param voxels: A 4-dimensional array of individual voxel dimensions in cMpc.
        :param cosmo: The EoRCosmology object for extracting the hubble constant.

        :returns (power, [kx, ky, kz]): The 3-dimensional power spectra of the data in The 3-dimensional power spectra of the data in mK^2 h^3 Mpc^-3 and the corresponding k-modes in three dimensions.
        """

        # apply the Blackman-Harris window in the frequency domain
        tdat = np.copy(data)
        for i in range(data.shape[2]):
            # Add half to centre the window
            tdat[:, :, i] = omath.blackman_harris_window(n=i + 0.5, N=data.shape[2]) * data[:, :, i]

        # Calculate temperature pertubation
        fperb = (tdat - tdat.mean()) * 1000  # mK

        # Calculate the total volume and voxel volume
        vol_vox = np.prod(voxels, axis=-1)
        vox_dim = np.mean(voxels, axis=(0, 1, 2), keepdims=False)

        # Fourier transform the pertubation on all dimensions
        ftdata = np.fft.fftn(fperb, norm="ortho")  # mK cMpc^-3

        # Take the square magnitude to remove complex values, multiply by voxel volume
        power = cosmo.little_h * vol_vox * np.abs(ftdata) ** 2  # mK^2 h^3 Mpc^-3

        # ki is different depending on dimension. Multiply by 2 pi: f -> k; multiply by little h; 1/cMpc -> h/Mpc
        ki = [
            2 * np.pi * np.fft.fftfreq(nx, d=vx) * cosmo.little_h for (nx, vx) in zip(list(data.shape), list(vox_dim))
        ]

        kx, ky, kz = np.meshgrid(*ki, indexing="ij")

        return power, [kx, ky, kz]

    @staticmethod
    def cylindrical_power_spectrum(power, ki, bins=[64, 32]):
        """Calculate the cylindrical power spectra for a given volume of data. Assumes that the voxel data is uniform.

        :param power: The 3-dimensional power spectra of the data in mK^2 h^3 Mpc^-3.
        :param ki: The corresponding [kx, ky, kz] values.
        :param bins: The number of k_perp and k_para bins as a tuple.

        :returns (power_2d, modal_agg, k_perp_edges, k_para_edges): The 2-dimensional power spectra, the edges of all of the k modes in the perpendicular and parallel dimensions.
        """

        # Extract ki
        [kx, ky, kz] = ki

        # Get the k perp and para modes
        k_perp = np.hypot(kx, ky)
        k_para = np.abs(kz)  # Remove negative values from kz

        # Calculate k-range max values using the 2 cMpc cell size
        k_para_max = k_para.max()  # cosmology.little_h * np.pi / 2
        k_perp_max = k_perp.max()  # np.sqrt(2) * cosmology.little_h * np.pi / 2

        # Calculate the minimum bin size and multiply by 2 pi (for the perp mode)
        k_para_min = k_para[np.nonzero(k_para)].min()  # cosmology.little_h / 128
        k_perp_min = k_perp[np.nonzero(k_perp)].min()  # np.sqrt(2) * cosmology.little_h / 1024

        # Collapse the power spectrum data into a 2D array
        power_agg, k_perp_edges, k_para_edges = np.histogram2d(
            k_perp.ravel(),
            k_para.ravel(),
            bins=bins,
            weights=power.ravel(),
            range=[[k_perp_min, k_perp_max], [k_para_min, k_para_max]],
        )

        # Get the number of modes used to aggregate in a specific bin
        modal_agg, _, _ = np.histogram2d(
            k_perp.ravel(),
            k_para.ravel(),
            bins=bins,
        )

        power_2d = np.divide(power_agg, modal_agg)

        return power_2d, ostr.mean_bins(k_perp_edges), ostr.mean_bins(k_para_edges)

    @staticmethod
    def dimensionless_power_spectrum(power, k_perp, k_para, bins=32, cosmo=None, z_centre=5.5):
        """Calculate the dimensionless, linear power spectra for a given 2DPS. Assumes that the voxel data is uniform.

        :param power: The 2-dimensional power spectra in mK^2 h^3 Mpc^-3.
        :param k_perp: The bin centres of the k modes in the perpendicular dimension.
        :param k_para: The bin centres of the k modes in the parallel dimension.
        :param bins: The number of k bins to use.
        :param cosmo: The EoRCosmology object for calculating the foreground wedge. If None, then no wedge is applied.
        :param z_centre: The central redshift value of the chunk.

        :returns (power_1d, k_edges): The dimensionless power spectra in mK^2, the edges of the k bins.
        """

        # Extract k-grid
        kpe, kpa = np.meshgrid(k_perp, k_para, indexing="ij")

        # Apply the foreground wedge
        if cosmo is not None:
            fw, _ = omath.foreground_wedge(kpe, kpa, cosmo=cosmo, z_centre=z_centre)
            pdwg = power * fw
        else:
            pdwg = power

        # Get the k modes
        k = np.hypot(kpe, kpa)

        # Calculate k-range max values using the 2 cMpc cell size
        k_max = k.max()  # np.sqrt(2) * cosmology.little_h * np.pi / 2

        # Calculate the minimum bin size and multiply by 2 pi (for the perp mode)
        k_min = k[np.nonzero(k)].min()  # np.sqrt(2) * cosmology.little_h / 1024

        # Get the dimensionless power
        pdim = pdwg * k**3 / (2 * np.pi**2)

        # Collapse the power spectrum data into a 2D array
        power_agg, k_edges = np.histogram(k.ravel(), bins=bins, weights=pdim.ravel(), range=[k_min, k_max])

        # Get the number of modes used to aggregate in a specific bin
        modal_agg, _ = np.histogram(k.ravel(), bins=bins, range=[k_min, k_max])

        power_1d = np.divide(power_agg, modal_agg)

        return power_1d, ostr.mean_bins(k_edges)

    def eval_power_spectrum(chunk, cvox):
        """Evaluate the power spectrum for a single redshift chunk.

        :param chunk: The array of chunk brightness temperature values.
        :param cvox: The array of chunk voxel dimensions.

        :return (ps2d, k_perp, k_para, ps1d, k): The 2D power spectra with its k_perp and k_para bins, as well as the 1D dimensionless power spectra with its k bins.
        """

        ps2d, k_perp, k_para = PSTools.cylindrical_power_spectrum(*PSTools.cubic_power_spectra(chunk, cvox))

        ps1d, k = PSTools.dimensionless_power_spectrum(ps2d, k_perp, k_para)

        return ps2d, k_perp, k_para, ps1d, k


class PSAnalysisPipeline:
    """Similar to the BTAnalysisPipeline but for Analysing CHIPS data, combines all of the power spectra analysis subtools."""
