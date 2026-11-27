#!.venv/bin/python

# ===

from matplotlib import pyplot as plt, colors as plc
import numpy as np
import sys
import h5py
from astropy import units as u, constants as c

# ===

# sys.path.append("/home/ohrf/Desktop/Curtin/phd-programs/")
sys.path.append("/scratch/mwaeor/ohrw/phd-programs")
from oskareor.reformatter import SimulationReformatter as simref
from oskareor.power_spectra_analysis import ChipsManager as chipsm, PSTools as pst
from chips_wrappers.ps_methods.chips_data import ChipsDataProducts as ChipsDataProducts
from chips_wrappers.setup.make_args import get_args, FakeArgs
from oskareor.skalow_calc import FileManager as ofmg, SKAMath as omath

# ===

path = ofmg.expand_path("/scratch/mwaeor/ohrw/oskareor.data/simulations/project1/fiducial_lightcone_simulation.h5")

data, z_ref, v, voxels, cosmology = simref.extract_h5_data(path)

cosmology.h0 = 67.66 * u.km / u.s / u.Mpc

# ===

# Import the h5 File Data
file = h5py.File(path)
data = np.array(file["lightcones"]["brightness_temp"])
dist = np.array(file["lightcone_distances"])
inpd = file["InputParameters"]

# ===

# Cosmology data
cosmo = cosmology
btzs = cosmology.dz_to_z(dist * u.Mpc)
print(btzs)

frqs = 1.420e9 / (btzs + 1)
print(frqs)

strs = np.arctan(2 / dist)
print(strs)

# ===

# Import the FITS data

from astropy.io import fits

print(data.shape)

with fits.open("/scratch/mwaeor/ohrw/oskareor.data/oskar.data/fiducial/fiducial_datacube.fits") as hdul:
    data = np.flip(np.moveaxis(hdul[0].data, 0, -1), -1)

print(data.shape)

# ===

# Calculate BT's from stokes I data provided by fits file
for i, (f, s) in enumerate(zip(frqs, strs)):
    data[:, :, i] = (omath.flux_to_brightness_temperature(fv=data[:, :, i], fxy=f, dphi=s, dtheta=s) * 1000) - 27

print(data.max(), data.min())

# ===

temps_unscaled = data.flatten()

zero_point = temps_unscaled.max()  # Define the zero point to force logarithmic bins
temps = zero_point - temps_unscaled
temps = temps[temps > 0]

# ===

# Print channel centres, widths, and start/end points and append to an accumulator list


def chunk_array(data, voxels, n_chunks=32, show_chunks=False, z_ref=5.5, cosmo=cosmology):
    """Split a large 3-dimensional datacube into a series of overlapping, evenly distributed, and fixed-size chunks for analysis in a power spectra.

    :param data: The 3-dimensional array to divide into chunks.
    :param voxels: The 4-dimensional voxel array to divide into chunks.
    :param n_chunks: The number of chunks to split the data in.
    :param show_chunks: Whether to show the chunk boundaries, defaults to False
    :param z_ref: The starting redshift of the box.
    :param cosmo: The EoRCosmology object for the box.

    :return (chunks, voxels, z_bins): Two 4-dimensional arrays containing each 3-dimensional chunk, one containing temperature values in Kelvin and the other voxel dimensions in cMpc. The third argument describes the start, middle, and end z-values for each chunk.
    """

    print("Dividing array into chunks ...")

    b = 2 * int(data.shape[2] / n_chunks)  # Chunk size
    dz_inc = np.mean(voxels[:, :, :, 2]) * u.Mpc  # Average voxel depth
    dz_ref = cosmo.z_to_dz(z_ref) + dz_inc

    chunk_values = []
    chunk_voxels = []
    chunk_z_bins = []

    if show_chunks:
        print("Chunk data:")
        print("+" + "—" * 86 + "+")
        print("| Chunk No.        Range ([start, end])    Centre      Length    Redshift Bins         |")
        print("+" + "—" * 86 + "+")

    for k in range(n_chunks - 1):
        c = int(b / 2 * (k + 1))
        i = int(c - b / 2)
        f = int(c + b / 2)

        z = np.array(
            [
                cosmo.dz_to_z(dz_ref + i * dz_inc),
                cosmo.dz_to_z(dz_ref + c * dz_inc),
                cosmo.dz_to_z(dz_ref + f * dz_inc),
            ]
        )

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
    chunk_z_bins = np.array(chunk_z_bins)

    if show_chunks:
        print("+" + "—" * 86 + "+")

    print(f"Final post-chunk shape: {chunk_values.shape}")

    return chunk_values, chunk_voxels, chunk_z_bins


chunk_values, chunk_voxels, chunk_z_bins = chunk_array(data=data, voxels=voxels, show_chunks=True, z_ref=z_ref)

# ===


def mean_bins(arr):
    """Calculate the arithmetic mean between each point in an array."""
    return (arr[1:] + arr[:-1]) / 2


def blackman_harris_window(n, N):
    """Three-term Blackman-Harris Windowing function."""
    a0 = 0.4243801
    a1 = 0.4973406
    a2 = 0.0782793

    return a0 - a1 * np.cos(2 * np.pi * n / N) + a2 * np.cos(4 * np.pi * n / N)


def foreground_wedge_gradient(cosmo, z_centre=5.5):
    """Evaluates the gradient of the foreground wedge."""
    dz = cosmo.z_to_dz(z_centre)
    fov = cosmo.fov_at_z(z_centre)
    Hz = cosmo.cosmo.efunc(z_centre) * cosmo.h0

    mz = (np.sin(fov) * dz * Hz / (c.c * (1 + z_centre))).to_value(u.dimensionless_unscaled)

    return mz


def foreground_wedge(k_perp, k_para, cosmo, z_centre=5.5):
    """Filter for discriminating out k-modes in the foreground. Assumes units of cMpc."""

    mz = foreground_wedge_gradient(cosmo=cosmo, z_centre=z_centre)
    k_compare = mz * k_perp

    return (k_para >= k_compare).astype(int), mz


# ===


def plot_2dps_axis(ax, k_perp, k_para, ps2d, z_bin, cosmo=None, plot_extras=True):
    img = ax.pcolormesh(k_perp, k_para, ps2d.T, cmap="viridis", norm="log")

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.tick_params(which="minor", length=3)
    ax.set_title(r"$z = " + f"{z_bin[1]:.2f}" + r"\pm" + f"{(z_bin[2] - z_bin[0])/2:.2f}" + r"$")

    if cosmo is not None:
        mz = foreground_wedge_gradient(cosmo=cosmology, z_centre=z_bin[1])
        plot_line_on_image(ax, x=k_perp, m=mz, b=0, label="Wedge Boundary")

    if plot_extras:
        ax.set_title(
            r"Pre-OSKAR Two-Dimensional Power Spectra for "
            + r"$z = "
            + f"{z_bin[1]:.2f}"
            + r"\pm"
            + f"{(z_bin[2] - z_bin[0])/2:.2f}"
            + r"$"
        )
        ax.set_xlabel(r"$k_{\perp}$ ($\mathrm{h\,Mpc}^{-1}$)")
        ax.set_ylabel(r"$k_{\parallel}$ ($\mathrm{h\,Mpc}^{-1}$)")

    return img


def plot_line_on_image(ax, x=None, m=0, b=0, label="", wedge=True):
    """Plots a line with formula y=mx+b ontop of an image."""

    xlim, ylim = ax.get_xlim(), ax.get_ylim()

    if x is None:
        x = np.linspace(0, 100, 100)
    x = np.asarray(x)

    y = m * x + b

    ax.plot(x, y, "k-", label=label)

    # If wedge is true shade the region below the line
    if wedge:
        ax.fill_between(x, ylim[0], y, color="black", alpha=0.4, label=label)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def plot_1dps_axis(ax, k, ps1d, z_bin, color=(0.0, 0.0, 0.0)):
    """Plots a 1-dimensional power spectrum onto a given axis."""

    ax.plot(
        mean_bins(k),
        ps1d,
        label=r"$z = " + f"{z_bin[1]:.2f}" + r"\pm" + f"{(z_bin[2] - z_bin[0])/2:.2f}" + r"$",
        color=color,
    )


# ===

fig, ax = plt.subplots(nrows=1, ncols=2, layout="constrained")

[ax1, ax2] = ax

tol_color = [
    "#332288",
    "#117733",
    "#44aa99",
    "#88ccee",
    "#ddcc77",
    "#cc6677",
    "#aa4499",
    "#882255",
]

for i in range(8):
    # Both
    chunk_sel = int(i * chunk_values.shape[0] / 8) + 2
    czbn = chunk_z_bins[chunk_sel]

    # Pre
    chunk = chunk_values[chunk_sel, :, :, :]
    cvox = chunk_voxels[chunk_sel, :, :, :]

    ps2d, k_perp, k_para = pst.cylindrical_power_spectrum(*pst.cubic_power_spectra(chunk, cvox))
    ps1d, k = pst.dimensionless_power_spectrum(ps2d, k_perp, k_para)

    plot_1dps_axis(ax1, k, ps1d, czbn, tol_color[i])

    # Post
    ps2d, k_perp, k_para = chipsm.extract_and_load_2dps(
        f"/scratch/mwaeor/ohrw/oskareor.data/oskar.data/fiducial/ps2ds/fiducial_ps2d_data_chunk_{chunk_sel}.npz"
    )
    # ps2d, k_perp, k_para = chipsm.extract_and_load_2dps(f"/scratch/mwaeor/ohrw/oskareor.data/2dps_npz_templates/point_ps2ds/point_ps2d_data_chunk_{chunk_sel}.npz")
    ps1d, k = pst.dimensionless_power_spectrum(ps2d, k_perp, k_para)

    plot_1dps_axis(ax2, k, ps1d, czbn, tol_color[i])


fig.supxlabel(r"$|\mathbf{k}|$ ($\mathrm{h} \mathrm{Mpc}^{-1}$)")
fig.supylabel(r"$\Delta_{21}^2 (k)$ ($\mathrm{mK}^2$)")

ax1.set_xscale("log")
ax2.set_xscale("log")
ax1.set_yscale("log")
ax2.set_yscale("log")

ax1.tick_params(which="minor", length=3)
ax2.tick_params(which="minor", length=3)

fig.legend(loc="center left", bbox_to_anchor=(1, 0.5))

fig.set_size_inches((12, 5))

ax1.set_title(r"Pre-OSKAR 1D Dimensionless Power Spectrum")
ax2.set_title(r"Post-OSKAR 1D Dimensionless Power Spectrum")

plt.savefig("img/chips_1dps.png")

# ===

fig = plt.figure(figsize=(12, 12), layout="constrained")
grid = fig.add_gridspec(5, 4, height_ratios=[1, 1, 0.12, 1, 1])
ax = np.array([fig.add_subplot(grid[row, col]) for row in (0, 1, 3, 4) for col in range(4)])
title_ax = fig.add_subplot(grid[2, :])
title_ax.set_axis_off()

# Pre-OSKAR
for i in range(8):
    chunk_sel = int(i * chunk_values.shape[0] / 8) + 2

    chunk = chunk_values[chunk_sel, :, :, :]
    cvox = chunk_voxels[chunk_sel, :, :, :]
    czbn = chunk_z_bins[chunk_sel]

    ps2d, k_perp, k_para = pst.cylindrical_power_spectrum(*pst.cubic_power_spectra(chunk, cvox))

    img = plot_2dps_axis(ax[i], k_perp, k_para, ps2d, czbn, cosmology, False)

# Post-OSKAR
for i in range(8):
    chunk_sel = int(i * chunk_values.shape[0] / 8) + 2
    czbn = chunk_z_bins[chunk_sel]

    ps2d, k_perp, k_para = chipsm.extract_and_load_2dps(
        f"/scratch/mwaeor/ohrw/oskareor.data/oskar.data/fiducial/ps2ds/fiducial_ps2d_data_chunk_{chunk_sel}.npz"
    )
    # ps2d, k_perp, k_para = chipsm.extract_and_load_2dps(f"/scratch/mwaeor/ohrw/oskareor.data/2dps_npz_templates/point_ps2ds/point_ps2d_data_chunk_{chunk_sel}.npz")

    img = plot_2dps_axis(ax[i + 8], k_perp, k_para, ps2d, czbn, cosmology, False)

# Add colorbar common to both sections
cb = fig.colorbar(
    img, ax=ax.tolist(), label=r"$P_{21}(k)$ ($\mathrm{mK}^2 \mathrm{h}^3 \mathrm{Mpc}^{-3}$)", extend="both"
)

title = fig.suptitle(r"Pre-OSKAR Two-Dimensional Power Spectra")
title_ax.text(
    0.5,
    0.5,
    r"Post-OSKAR Two-Dimensional Power Spectra",
    ha="center",
    va="center",
    fontproperties=title.get_fontproperties(),
)
fig.supxlabel(r"$k_{\perp}$ ($\mathrm{h\,Mpc}^{-1}$)")
fig.supylabel(r"$k_{\parallel}$ ($\mathrm{h\,Mpc}^{-1}$)")

plt.savefig("img/chips_2dps.png")
