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

            chunk_bounds.append(np.array([i, c, f]))
            chunk_z_bins.append(z)

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
                15.0,
                c.c.si.value,
            ],
            check=True,
        )

        subprocess.run(
            [
                "/software/projects/mwaeor/ctrott/setonix/chips_2025_ska/prepare_ska",
                "side",
                n_channels,
                0,
                "'yy'",
                "side",
                1,
            ],
            check=True,
        )

        subprocess.run(
            [
                "/software/projects/mwaeor/ctrott/setonix/chips_2025_ska/lssa_fg_ska",
                "side",
                n_channels,
                k_perp_bins,
                "'yy'",
                u_max,
                "chips.out",
                0,
                1,
            ],
            check=True,
        )

    @staticmethod
    def extract_and_save_2dps(outdir, temp_parent_dir="."):
        """Extracts the CHIPS output of a chunk, combines the data, and then saves it to a .npz file for later useage."""

        # FIXME: Move chips output to safe location

    @staticmethod
    def run_chips_on_model(
        uvfits_file,
        temp_parent_dir=".",
        outdir=".",
        voxel_depth=2,
        n_channels=64,
        n_chunks=32,
        k_perp_bins=64,
        u_max=500.0,
    ):
        """Runs the entirety of the CHIPS model."""

        print("Creating the output directory ./PSA ...")

        # Create the directory and/or clear it
        subprocess.run(["rm", "-rf", temp_parent_dir + "/PSA/chips_output"], check=True)
        subprocess.run(["mkdir", "-p", temp_parent_dir + "/PSA/chips_output"], check=True)

        print("Generating chunks ...")

        # Create an array of chunks
        chunk_bins = Chunker.generate_chunks(
            n_channels=n_channels, voxel_depth=voxel_depth, n_chunks=n_chunks, show_chunks=False
        )

        # Create the temporary UVFits file location for CHIPS to analyse
        uvfits_loc = temp_parent_dir + ChipsManager.DEFAULT_TEMP_FILE

        for i, chunk_bounds in enumerate(chunk_bins):
            print(f"Extracting and evaluating chunk: {i} of {n_chunks}")

            ChipsManager.extract_uvfits_chunk(
                uvfits_file=uvfits_file, chunk_boundaries=chunk_bounds, outfile=uvfits_loc
            )

            ChipsManager.run_chips(
                uvfits_loc,
                temp_parent_dir=temp_parent_dir + "/PSA/chips_output/",
                n_channels=n_channels,
                k_perp_bins=k_perp_bins,
                u_max=u_max,
            )

            ChipsManager.extract_and_save_2dps(outdir, temp_parent_dir + "/PSA/chips_output")


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

        return power_2d, k_perp_edges, k_para_edges

    @staticmethod
    def dimensionless_power_spectrum(power, k_perp_edges, k_para_edges, bins=32, cosmo=None, z_centre=5.5):
        """Calculate the dimensionless, linear power spectra for a given 2DPS. Assumes that the voxel data is uniform.

        :param power: The 2-dimensional power spectra in mK^2 h^3 Mpc^-3.
        :param k_perp_edges: The edges of all of the k modes in the perpendicular dimension.
        :param k_perp_edges: The edges of all of the k modes in the parallel dimension.
        :param bins: The number of k bins to use.
        :param cosmo: The EoRCosmology object for calculating the foreground wedge. If None, then no wedge is applied.
        :param z_centre: The central redshift value of the chunk.

        :returns (power_1d, k_edges): The dimensionless power spectra, the edges of the k bins.
        """

        # Extract k-grid
        k_perp, k_para = np.meshgrid(ostr.mean_bins(k_perp_edges), ostr.mean_bins(k_para_edges), indexing="ij")

        # Apply the foreground wedge
        if cosmo is not None:
            fw, _ = omath.foreground_wedge(k_perp, k_para, cosmo=cosmo, z_centre=z_centre)
            pdwg = power * fw
        else:
            pdwg = power

        # Get the k modes
        k = np.hypot(k_perp, k_para)

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

        return power_1d, k_edges


class PSAnalysisPipeline:
    """Similar to the BTAnalysisPipeline but for Analysing CHIPS data, combines all of the power spectra analysis subtools."""
