#!.venv/bin/python

# ===

# Import python modules
import os
import sys
import numpy as np
from matplotlib import pyplot as plt, colors
from astropy import units as u, constants as c
from astropy.cosmology import FlatLambdaCDM as fmodel, z_at_value as getz
import h5py

# import py21cmfast as p21c
import json

# Import local oskareor models
sys.path.append("/scratch/mwaeor/ohrw/phd-programs")
from oskareor.skalow_calc import SKAMath as omath, EoRCosmology as Cosmo, FileManager as ofmg

# ===

# Import the h5 File Data
path = ofmg.expand_path("/scratch/mwaeor/ohrw/oskareor.data/simulations/project1/fiducial_lightcone_simulation.h5")
file = h5py.File(path)
data = np.array(file["lightcones"]["brightness_temp"])
dist = np.array(file["lightcone_distances"])
inpd = file["InputParameters"]

# ===

print(list(file))
print(dist)
print(list(file["lightcones"]))

print(file["InputParameters"]["simulation_options"].attrs["HII_DIM"])

# ===


def hdf5_to_dict(grp):
    """Load all data from an HDF5 Group into a dict.

    Essentially the same as toml.load() but for HDF5.
    """
    out = dict(grp.attrs)

    for k, v in grp.items():
        if isinstance(v, h5py.Group):
            out[k] = hdf5_to_dict(v)
        else:
            out[k] = v[()]

    return out


print(json.dumps(hdf5_to_dict(inpd), indent=4, sort_keys=True, default=str))
inpd = hdf5_to_dict(inpd)

# ===

# Cosmology data
cosmo = Cosmo(
    inpd["cosmo_params"]["hlittle"] * 100 * u.km / u.s / u.Mpc, inpd["cosmo_params"]["OMm"], inpd["cosmo_params"]["OMb"]
).cosmo
btzs = getz(cosmo.comoving_distance, dist * u.Mpc).to_value()
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

zero_point = temps.max()  # Define the zero point to force logarithmic bins
temps = zero_point - temps_unscaled
temps = temps[temps > 0]

# ===

fig, ax = plt.subplots(layout="constrained")

ax.set_title("Full-Cube Temperature Distribution")
ax.set_xlabel(r"$T_b$ (mK)")
ax.set_ylabel(r"Frequency")

bins = np.logspace(np.log10(temps.min()), np.log10(temps.max()), 50)
ax.hist(temps, bins=bins, histtype="bar")

ax.xaxis.set_inverted(True)
ax.set_xscale("log")
ax.set_yscale("log")

print(np.trunc(-np.log10(-temps_unscaled.min())), np.trunc(np.log10(temps_unscaled.max())))

major_tick_base = [
    np.arange(3.0, 0.0, -1.0),
    0,
    np.arange(1.0, 2.0, 1.0),
]

print(major_tick_base)

major_ticks = np.concatenate(
    [
        zero_point + 10 ** major_tick_base[0],
        [zero_point],
        zero_point - 10 ** major_tick_base[2],
    ]
)

major_tick_labels = np.concatenate(
    [
        np.array([r"$-10^" + str(int(x)) + r"$" for x in major_tick_base[0]]),
        ["0"],
        np.array([r"$10^" + str(int(x)) + r"$" for x in major_tick_base[2]]),
    ]
)

print(major_ticks, major_tick_labels)

ax.set_xticks(major_ticks)
ax.set_xticklabels(major_tick_labels)

ax.tick_params(which="minor", length=3)
# ax.xaxis.set_major_locator(plt.LogLocator(base=10))

plt.savefig("img/temps.png")

# ===

# EoR colour map
EoR_colour = colors.LinearSegmentedColormap.from_list(
    "mycmap",
    [(0, "white"), (0.125, "yellow"), (0.25, "orange"), (0.325, "red"), (0.5, "black"), (0.75, "blue"), (1, "cyan")],
)

eor_cmap = {
    "map": EoR_colour,
    "norm": colors.TwoSlopeNorm(0, -160, 40),
    "ticks": np.arange(-160, 41, 10),
    "bounds": np.linspace(-160, 40, 257),
    "dist": "uniform",
}

# Snapshots

channels = 126 + np.arange(8) * 128

fig = plt.figure(layout="constrained", figsize=(12, 12))
gs = fig.add_gridspec(nrows=4, ncols=4)
ax = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(4)]

fig.suptitle("Single-Channel Slices for the Fiducial Simulation Cube")

for i, chn in enumerate(channels):
    im = ax[i].imshow(data[:, :, chn], cmap="inferno")
    ax[i].set_xticks([])
    ax[i].set_yticks([])
    ax[i].set_title(f"Ch. #{chn+1} (z = {btzs[chn]:3.2f})")

# Lightcone


def custom_image_ticks(values, labels, scaling=None):
    if scaling is None:
        scaling = lambda x: x
    scaling = np.vectorize(scaling)
    values = scaling(values)
    labels = scaling(labels)
    scales = np.interp(labels, values, np.arange(values.shape[0]))
    return scales, labels


axl = fig.add_subplot(gs[2:, :])

axl.set_title(r"Lightcone of the Fiducial Model at $0^{\circ}\,\mathrm{RA}$")

im = axl.imshow(data[255, :, :], cmap="inferno")

axl.set_xticks(*custom_image_ticks(btzs, np.arange(6, 14, 1)))
axl.set_yticks(*custom_image_ticks(np.arange(512)*2, np.arange(0, 1025, 128), scaling=(lambda x: x - 512)))

nrows, ncols = im.get_size()
aspect_ratio = ncols / nrows

axl.set_ylabel(r"Physical Distance from Zenith $\sim -27^{\circ}\,\mathrm{Dec}$ ($\mathrm{cMpc}$)")
axl.set_xlabel("Redshift")

cbar = fig.colorbar(
    im, ax=[*ax, axl], label="Brightness Temperature (mK)", aspect=40
)  # , boundaries=eor_cmap["bounds"], spacing=eor_cmap["dist"], ticks=eor_cmap["ticks"])

cbar.ax.yaxis.set_ticks_position("left")
cbar.ax.yaxis.set_label_position("left")

plt.savefig("img/bta_fits_cube.png")
