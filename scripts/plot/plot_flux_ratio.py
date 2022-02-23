import os
os.chdir("/home/brian/Degrees/masters/code/python/projects/flux_suppression")

import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from astropy.io import fits
import astLib.astWCS
import Tigger
from tqdm import tqdm


# Metric Function - Root Mean Square Error
def rmse(A, B):
    res = A - B
    return np.sqrt(np.vdot(res, res).real/res.size)

# Params for plots
rcParams.update({'font.size': 11, 'font.family': 'serif'})

arrays = {
      "kat7": {
            1: [
                  (100, -6, -1)
            ],

            2: [
                  (100, -6, -2),
                  (70, -6, -2.7)
            ]

      },
      "vlab": {
            1: [
                  (100, -6, -1)
            ],

            2: [
                  (100, -5, -1.7),
                  (70, -6, -3)
            ],

            100: [
                  (100, -6, -1),
                  (70, -6, -1),
                  (50, -6, -2),
                  (30, -6, -3)
            ]
      },
      "meerkat": {
            1: [
                  (100, -6, -1)
            ],

            2: [
                  (100, -5, -1.7),
                  (70, -5, -3.7)
            ],

            100: [
                  (100, -6, -1),
                  (70, -6, -1),
                  (50, -6, -2.5),
                  (30, -6, -3.2)
            ]
      }
}

tina_path = "/net/tina/vault2-tina/welman/flux_suppression/"

# Options
MODE = "image"
FORMAT = ".fits"
n_points = 16

# All time solution intervals
tints = np.arange(12, 252, 12)
tints = np.insert(tints, 0, 8)
tints = np.insert(tints, 0, 4)
tints = np.insert(tints, 0, 1)

kalcal_flux = np.zeros((3, 3, 4, n_points), dtype=np.float64)
kalcal_sources = np.zeros((2, 4, 2, 100), dtype=np.float64)
kalcal_sources = {
    "vlab": {
        100: {
            "flux": None,
            "true_flux": None
        },

        70: {
            "flux": None,
            "true_flux": None
        },

        50: {
            "flux": None,
            "true_flux": None
        },

        30: {
            "flux": None,
            "true_flux": None
        }
    },

    "meerkat": {
        100: {
            "flux": None,
            "true_flux": None
        },

        70: {
            "flux": None,
            "true_flux": None
        },

        50: {
            "flux": None,
            "true_flux": None
        },

        30: {
            "flux": None,
            "true_flux": None
        }
    }
}
with tqdm(total=272) as pbar:
      for i, (ant, sources) in enumerate(arrays.items()):
        for j, (source, setup) in enumerate(sources.items()):
            for k, (percent, lb, ub) in enumerate(setup):
                sigma_fs = np.round(np.logspace(lb, ub, n_points), 7)
                ADJ = 0
                if source == 1 or source == 2:
                    if source == 2 and int(percent) == 70:
                        ADJ = 1
                    SKY_SUFFIX = f"-src-{source - ADJ}.txt"
                else:
                    SKY_SUFFIX = f"-fp-{percent}-src-{source}.txt"

                true_flux = [src.flux.I for src in Tigger.load(f"skymodels/{ant}/model" + SKY_SUFFIX, verbose=False)]

                for l, sigma_f in enumerate(sigma_fs):
                    restored_image = tina_path + f"images/kalcal-img/{ant}/src-{source}/fp-{percent}/sigma_{sigma_f}/{ant}-fp-{percent}-src-{source}-{sigma_f}-{MODE + FORMAT}"
                    flux = get_flux(f"skymodels/{ant}/model" + SKY_SUFFIX, restored_image, 5)

                    if source == 100:
                        kalcal_flux[i, j, k, l] = np.mean(flux[:5]/true_flux[:5])
                        kalcal_sources[ant][percent]["true_flux"] = true_flux
                        kalcal_sources[ant][percent]["flux"] = flux
                    else:
                        kalcal_flux[i, j, k, l] = np.mean(flux/true_flux)
                    pbar.update(1)

# Plot RMSE
print("# Plotting figure 1")
fig, axes = plt.subplots(1, 3, figsize=(16, 9))

# Labels
axes[0].set_ylabel("Flux Ratio")
axes[0].set_xlabel(r"Process Noise - $\sigma_f$")
axes[1].set_xlabel(r"Process Noise - $\sigma_f$")
axes[2].set_xlabel(r"Process Noise - $\sigma_f$")

axes[0].set_title(f"KAT-7")
axes[1].set_title(f"VLA-B")
axes[2].set_title(f"MeerKAT")

fig.suptitle(f"Calibration over Single Source\nFlux Ratio between True Image and Restored Image (Kalcal over Process Noise)")

ones = np.ones(n_points, dtype=np.float64)

# Plot values
sigma_fs = np.round(np.logspace(arrays["kat7"][1][0][1], arrays["kat7"][1][0][2], n_points), 7)
axes[0].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[0].plot(sigma_fs, kalcal_flux[0, 0, 0], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[0].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["vlab"][1][0][1], arrays["vlab"][1][0][2], n_points), 7)
axes[1].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[1].plot(sigma_fs, kalcal_flux[1, 0, 0], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[1].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["meerkat"][1][0][1], arrays["meerkat"][1][0][2], n_points), 7)
axes[2].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[2].plot(sigma_fs, kalcal_flux[2, 0, 0], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[2].set_xscale("log")

# Set legend
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, facecolor='white', loc='lower center', ncol=2)

# Change x-ticks
for ax in axes:
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.4g'))

# Save plot
print("# Saving figure 1")
plt.savefig(f"plots/comparison/flux-ratio-all-src-1.png")

# Plot RMSE
print("# Plotting figure 2")
fig, axes = plt.subplots(2, 3, figsize=(16, 9))

# Labels
axes[0, 0].set_ylabel("Flux Ratio")
axes[1, 0].set_ylabel("Flux Ratio")
axes[1, 0].set_xlabel(r"Process Noise - $\sigma_f$")
axes[1, 1].set_xlabel(r"Process Noise - $\sigma_f$")
axes[1, 2].set_xlabel(r"Process Noise - $\sigma_f$")

axes[0, 0].set_title(f"KAT-7\n100% Modelled Flux")
axes[0, 1].set_title(f"VLA-B\n100% Modelled Flux")
axes[0, 2].set_title(f"MeerKAT\n100% Modelled Flux")
axes[1, 0].set_title(f"70% Modelled Flux")
axes[1, 1].set_title(f"70% Modelled Flux")
axes[1, 2].set_title(f"70% Modelled Flux")

fig.suptitle(f"Calibration over Two Sources\nFlux Ratio between True Image and Restored Image (Kalcal over Process Noise)")

# Plot values

sigma_fs = np.round(np.logspace(arrays["kat7"][2][0][1], arrays["kat7"][2][0][2], n_points), 7)
axes[0, 0].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[0, 0].plot(sigma_fs, kalcal_flux[0, 1, 0], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[0, 0].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["vlab"][2][0][1], arrays["vlab"][2][0][2], n_points), 7)
axes[0, 1].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[0, 1].plot(sigma_fs, kalcal_flux[1, 1, 0], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[0, 1].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["meerkat"][2][0][1], arrays["meerkat"][2][0][2], n_points), 7)
axes[0, 2].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[0, 2].plot(sigma_fs, kalcal_flux[2, 1, 0], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[0, 2].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["kat7"][2][1][1], arrays["kat7"][2][1][2], n_points), 7)
axes[1, 0].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[1, 0].plot(sigma_fs, kalcal_flux[0, 1, 1], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[1, 0].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["vlab"][2][1][1], arrays["vlab"][2][1][2], n_points), 7)
axes[1, 1].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[1, 1].plot(sigma_fs, kalcal_flux[1, 1, 1], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[1, 1].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["meerkat"][2][1][1], arrays["meerkat"][2][1][2], n_points), 7)
axes[1, 2].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[1, 2].plot(sigma_fs, kalcal_flux[2, 1, 1], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[1, 2].set_xscale("log")

# Set legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, facecolor='white', loc='lower center', ncol=2)

for ax1 in axes:
    for ax in ax1:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.4g'))

# Save plot
print("# Saving figure 2")
plt.savefig(f"plots/comparison/flux-ratio-all-src-2.png")

# Plot RMSE
print("# Plotting figure 3")
fig, axes = plt.subplots(2, 2, figsize=(16, 9))

# Labels
axes[0, 0].set_ylabel("Flux Ratio")
axes[1, 0].set_ylabel("Flux Ratio")
axes[1, 0].set_xlabel(r"Process Noise - $\sigma_f$")
axes[1, 1].set_xlabel(r"Process Noise - $\sigma_f$")

axes[0, 0].set_title(f"100% Modelled Flux")
axes[0, 1].set_title(f"69% Modelled Flux")
axes[1, 0].set_title(f"47% Modelled Flux")
axes[1, 1].set_title(f"28% Modelled Flux")

fig.suptitle(f"VLA-B calibration over 100 Sources\nFlux Ratio between True Image and Restored Image (Kalcal over Process Noise)")

# Plot values

sigma_fs = np.round(np.logspace(arrays["vlab"][100][0][1], arrays["vlab"][100][0][2], n_points), 7)
axes[0, 0].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[0, 0].plot(sigma_fs, kalcal_flux[1, 2, 0], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[0, 0].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["vlab"][100][1][1], arrays["vlab"][100][1][2], n_points), 7)
axes[0, 1].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[0, 1].plot(sigma_fs, kalcal_flux[1, 2, 1], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[0, 1].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["vlab"][100][2][1], arrays["vlab"][100][2][2], n_points), 7)
axes[1, 0].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[1, 0].plot(sigma_fs, kalcal_flux[1, 2, 2], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[1, 0].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["vlab"][100][3][1], arrays["vlab"][100][3][2], n_points), 7)
axes[1, 1].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[1, 1].plot(sigma_fs, kalcal_flux[1, 2, 3], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[1, 1].set_xscale("log")

# Set legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, facecolor='white', loc='lower center', ncol=2)

for ax1 in axes:
    for ax in ax1:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.4g'))

# Save plot
print("# Saving figure 3")
plt.savefig(f"plots/comparison/flux-ratio-vlab-src-100.png")

# Plot RMSE
print("# Plotting figure 4")
fig, axes = plt.subplots(2, 2, figsize=(16, 9))

# Labels
axes[0, 0].set_ylabel("Flux Ratio")
axes[1, 0].set_ylabel("Flux Ratio")
axes[1, 0].set_xlabel(r"Process Noise - $\sigma_f$")
axes[1, 1].set_xlabel(r"Process Noise - $\sigma_f$")

axes[0, 0].set_title(f"100% Modelled Flux")
axes[0, 1].set_title(f"69% Modelled Flux")
axes[1, 0].set_title(f"47% Modelled Flux")
axes[1, 1].set_title(f"28% Modelled Flux")

fig.suptitle(f"MeerKAT calibration over 100 Sources\nFlux Ratio between True Image and Restored Image (Kalcal over Process Noise)")

# Plot values
sigma_fs = np.round(np.logspace(arrays["meerkat"][100][0][1], arrays["meerkat"][100][0][2], n_points), 7)
axes[0, 0].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[0, 0].plot(sigma_fs, kalcal_flux[2, 2, 0], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[0, 0].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["meerkat"][100][1][1], arrays["meerkat"][100][1][2], n_points), 7)
axes[0, 1].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[0, 1].plot(sigma_fs, kalcal_flux[2, 2, 1], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[0, 1].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["meerkat"][100][2][1], arrays["meerkat"][100][2][2], n_points), 7)
axes[1, 0].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[1, 0].plot(sigma_fs, kalcal_flux[2, 2, 2], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[1, 0].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["meerkat"][100][3][1], arrays["meerkat"][100][3][2], n_points), 7)
axes[1, 1].plot(sigma_fs, ones, linestyle="-", marker="s", color="black", label="True")
axes[1, 1].plot(sigma_fs, kalcal_flux[2, 2, 3], linestyle="-", marker="x", color="forestgreen", label="Kalcal")
axes[1, 1].set_xscale("log")

# Set legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, facecolor='white', loc='lower center', ncol=2)

for ax1 in axes:
    for ax in ax1:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.4g'))

# Save plot
print("# Saving figure 4")
plt.savefig(f"plots/comparison/flux-ratio-meerkat-src-100.png")

# Plot RMSE
print("# Plotting figure 5")
fig, axes = plt.subplots(2, 2, figsize=(16, 9))

# Labels
axes[0, 0].set_ylabel("Flux (Jy)")
axes[1, 0].set_ylabel("Flux (Jy)")
axes[1, 0].set_xlabel("Source Index")
axes[1, 1].set_xlabel("Source Index")

axes[0, 0].set_title(f"100% Modelled Flux")
axes[0, 1].set_title(f"69% Modelled Flux")
axes[1, 0].set_title(f"47% Modelled Flux")
axes[1, 1].set_title(f"28% Modelled Flux")

fig.suptitle(f"VLA-B calibration over 100 Sources\nFlux Ratio between True Image and Restored Image (Kalcal over Process Noise)")

# Plot values
axes[0, 0].plot(np.arange(len(kalcal_sources["vlab"][100]["true_flux"])) + 1, kalcal_sources["vlab"][100]["true_flux"], "+", color="black", label="True Sources")
axes[0, 0].plot(np.arange(len(kalcal_sources["vlab"][100]["true_flux"])) + 1, kalcal_sources["vlab"][100]["flux"], "-", color="forestgreen", label="Kalcal Sources")

axes[0, 1].plot(np.arange(len(kalcal_sources["vlab"][70]["true_flux"])) + 1, kalcal_sources["vlab"][70]["true_flux"], "+", color="black", label="True Sources")
axes[0, 1].plot(np.arange(len(kalcal_sources["vlab"][70]["true_flux"])) + 1, kalcal_sources["vlab"][70]["flux"], "-", color="forestgreen", label="Kalcal Sources")

axes[1, 0].plot(np.arange(len(kalcal_sources["vlab"][50]["true_flux"])) + 1, kalcal_sources["vlab"][50]["true_flux"], "+", color="black", label="True Sources")
axes[1, 0].plot(np.arange(len(kalcal_sources["vlab"][50]["true_flux"])) + 1, kalcal_sources["vlab"][50]["flux"], "-", color="forestgreen", label="Kalcal Sources")

axes[1, 1].plot(np.arange(len(kalcal_sources["vlab"][30]["true_flux"])) + 1, kalcal_sources["vlab"][30]["true_flux"], "+", color="black", label="True Sources")
axes[1, 1].plot(np.arange(len(kalcal_sources["vlab"][30]["true_flux"])) + 1, kalcal_sources["vlab"][30]["flux"], "-", color="forestgreen", label="Kalcal Sources")

# Set legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, facecolor='white', loc='lower center', ncol=2)

for ax1 in axes:
    for ax in ax1:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.4g'))

# Save plot
print("# Saving figure 5")
plt.savefig(f"plots/comparison/source-compare-vlab-src-100.png")

# Plot RMSE
print("# Plotting figure 6")
fig, axes = plt.subplots(2, 2, figsize=(16, 9))

# Labels
axes[0, 0].set_ylabel("Flux (Jy)")
axes[1, 0].set_ylabel("Flux (Jy)")
axes[1, 0].set_xlabel("Source Index")
axes[1, 1].set_xlabel("Source Index")

axes[0, 0].set_title(f"100% Modelled Flux")
axes[0, 1].set_title(f"69% Modelled Flux")
axes[1, 0].set_title(f"47% Modelled Flux")
axes[1, 1].set_title(f"28% Modelled Flux")

fig.suptitle(f"MeerKAT calibration over 100 Sources\nFlux Ratio between True Image and Restored Image (Kalcal over Process Noise)")

# Plot values
axes[0, 0].plot(np.arange(len(kalcal_sources["meerkat"][100]["true_flux"])) + 1, kalcal_sources["meerkat"][100]["true_flux"], "+", color="black", label="True Sources")
axes[0, 0].plot(np.arange(len(kalcal_sources["meerkat"][100]["true_flux"])) + 1, kalcal_sources["meerkat"][100]["flux"], "-", color="forestgreen", label="Kalcal Sources")

axes[0, 1].plot(np.arange(len(kalcal_sources["meerkat"][70]["true_flux"])) + 1, kalcal_sources["meerkat"][70]["true_flux"], "+", color="black", label="True Sources")
axes[0, 1].plot(np.arange(len(kalcal_sources["meerkat"][70]["true_flux"])) + 1, kalcal_sources["meerkat"][70]["flux"], "-", color="forestgreen", label="Kalcal Sources")

axes[1, 0].plot(np.arange(len(kalcal_sources["meerkat"][50]["true_flux"])) + 1, kalcal_sources["meerkat"][50]["true_flux"], "+", color="black", label="True Sources")
axes[1, 0].plot(np.arange(len(kalcal_sources["meerkat"][50]["true_flux"])) + 1, kalcal_sources["meerkat"][50]["flux"], "-", color="forestgreen", label="Kalcal Sources")

axes[1, 1].plot(np.arange(len(kalcal_sources["meerkat"][30]["true_flux"])) + 1, kalcal_sources["meerkat"][30]["true_flux"], "+", color="black", label="True Sources")
axes[1, 1].plot(np.arange(len(kalcal_sources["meerkat"][30]["true_flux"])) + 1, kalcal_sources["meerkat"][30]["flux"], "-", color="forestgreen", label="Kalcal Sources")

# Set legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, facecolor='white', loc='lower center', ncol=2)

for ax1 in axes:
    for ax in ax1:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.4g'))

# Save plot
print("# Saving figure 6")
plt.savefig(f"plots/comparison/source-compare-meerkat-src-100.png")

# Show plots
# plt.show()
