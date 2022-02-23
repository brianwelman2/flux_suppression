from casacore.tables import table
import os
os.chdir("/home/welman/masters/projects/flux_suppression")

import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from tqdm import tqdm
from africanus.calibration.utils import corrupt_vis


# Params for plots
rcParams.update({'font.size': 12, 'font.family': 'serif'})

# Metric Function - Root Mean Square Error
def mse(A, B):
    res = A - B
    return np.vdot(res, res).real/res.size

# Give sigma_n
sigma_n = np.sqrt(2)
weight = np.sqrt(1.0/2 * sigma_n**2)
W = 1.0/(np.sqrt(2) * sigma_n)

# All sigma_fs
#sigma_fs = np.round(np.hstack((np.arange(0.0005, 0.001, 0.0001), np.arange(0.001, 0.02, 0.001), np.arange(0.02, 1.02, 0.02))), 4)

tina_path = "/net/tina/vault2-tina/welman/flux_suppression/"

arrays = {
      "kat7": {
            1: [
                  (100, -6, -1)
            ],

            2: [
                  (100, -2, 2),
                  (70, -2, 2)
            ]

      },
      "vlab": {
            1: [
                  (100, -9, -3)
            ],

            2: [
                  (100, -3, 2),
                  (70, -3, 2)
            ],

            100: [
                  (100, -3, 1),
                  (70, -3, 1),
                  (50, -3, 1),
                  (30, -3, 1)
            ]
      },
      "meerkat": {
            1: [
                  (100, -9, -3)
            ],

            2: [
                  (100, -3, 1),
                  (70, -3, 1)
            ],

            100: [
                  (100, -3, 1),
                  (70, -3, 1),
                  (50, -3, 1),
                  (30, -3, 1)
            ]
      }
}

# Minimums
n_points = 16
mins = np.zeros((3, 3, 4, 2), dtype=np.float64)
kalcal_rmse = np.zeros((3, 3, 4, n_points, 4), dtype=np.float64)
decimals = 16

with tqdm(total=272) as pbar:
      for i, (ant, sources) in enumerate(arrays.items()):
            with table(f"ms/{ant}.ms", ack=False) as tb:
                  _, tbin_indices, tbin_counts = np.unique(tb.getcol("TIME"), return_index=True, return_counts=True)
                  ant1 = tb.getcol("ANTENNA1")
                  ant2 = tb.getcol("ANTENNA2")
                  for j, (source, setup) in enumerate(sources.items()):
                        for k, (percent, lb, ub) in enumerate(setup):
                              sigma_fs = np.round(np.logspace(lb, ub, n_points), decimals)
                              true_gains = np.load(tina_path + f"gains/true/{ant}/{ant}-src-{source}-gains.npy")[..., 0]

                              n_time, n_ant, n_chan, n_dir = true_gains.shape
                              n_row = n_time * (n_ant * (n_ant - 1)//2)
                              ADJ = 0

                              if source == 1 or source == 2:
                                    MODEL_PREFIX = ""
                                    DATA_PREFIX = ""
                                    if source == 2 and int(percent) == 70:
                                          ADJ = 1
                              else:
                                    MODEL_PREFIX = f"FP{percent}_"
                                    DATA_PREFIX = f"FP100_"

                              # model = tb.getcol(MODEL_PREFIX + f"SRC{source - ADJ}_MODEL")[..., np.newaxis, (0, 3)]
                              # clean_vis = tb.getcol(DATA_PREFIX + f"SRC{source}_DATA")[..., (0, 3)]

                              for l, sigma_f in enumerate(sigma_fs):
                                    filter_path = tina_path + f"gains/kalcal/{ant}/src-{source}/fp-{percent}/sigma_{sigma_f}/full-filter.npy"
                                    filter_gains = np.load(filter_path)[..., 0, 0]
                                    kalcal_rmse[i, j, k, l, 0] = mse(true_gains, filter_gains)

                                    smoother_path = tina_path + f"gains/kalcal/{ant}/src-{source}/fp-{percent}/sigma_{sigma_f}/full-smoother.npy"
                                    smoother_gains = np.load(smoother_path)[..., 0, 0]
                                    kalcal_rmse[i, j, k, l, 1] = mse(true_gains, smoother_gains)

                                    # corr_vis = corrupt_vis(tbin_indices, tbin_counts, ant1, ant2, filter_gains, model)
                                    # kalcal_rmse[i, j, k, l, 2] = rcs(clean_vis, corr_vis, tbin_indices, tbin_counts,
                                    #                                     ant1, ant2, smoother_gains, n_time, n_chan, n_ant, n_row)

                                    # corr_vis = corrupt_vis(tbin_indices, tbin_counts, ant1, ant2, smoother_gains, model)
                                    # kalcal_rmse[i, j, k, l, 3] = rcs(clean_vis, corr_vis, tbin_indices, tbin_counts,
                                    #                                     ant1, ant2, smoother_gains, n_time, n_chan, n_ant, n_row)

                                    pbar.update(1)


                              mins[i, j, k, 0] = sigma_fs[kalcal_rmse[i, j, k, :, 1].argmin()]
                              mins[i, j, k, 1] = kalcal_rmse[i, j, k, :, 1].min()

# Plot RMSE
print("# Plotting figure 1")
fig, axes = plt.subplots(1, 3, figsize=(16, 9))

# Labels
axes[0].set_ylabel("RMSE")
axes[0].set_xlabel(r"Process Noise - $\sigma_f$")
axes[1].set_xlabel(r"Process Noise - $\sigma_f$")
axes[2].set_xlabel(r"Process Noise - $\sigma_f$")

axes[0].set_title(f"KAT-7")
axes[1].set_title(f"VLA-B")
axes[2].set_title(f"MeerKAT")

fig.suptitle(f"Calibration over Single Source\nRMSE between True Gains and Calibrated Gains (Kalcal over Process Noise)")

# # Plot values
sigma_fs = np.round(np.logspace(arrays["kat7"][1][0][1],
                                arrays["kat7"][1][0][2], n_points), decimals)
axes[0].plot(sigma_fs, kalcal_rmse[0, 0, 0, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[0].plot(sigma_fs, kalcal_rmse[0, 0, 0, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[0].plot(sigma_fs, kalcal_rmse[0, 0, 0, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[0].plot(sigma_fs, kalcal_rmse[0, 0, 0, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[0].plot(mins[0, 0, 0, 0], mins[0, 0, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[0].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["vlab"][1][0][1],
                                arrays["vlab"][1][0][2], n_points), decimals)
axes[1].plot(sigma_fs, kalcal_rmse[1, 0, 0, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[1].plot(sigma_fs, kalcal_rmse[1, 0, 0, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[1].plot(sigma_fs, kalcal_rmse[1, 0, 0, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[1].plot(sigma_fs, kalcal_rmse[1, 0, 0, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[1].plot(mins[1, 0, 0, 0], mins[1, 0, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[1].set_xscale("log")
print(mins[:, 0, 0])
exit()
sigma_fs = np.round(np.logspace(arrays["meerkat"][1][0][1],
                                arrays["meerkat"][1][0][2], n_points), decimals)
axes[2].plot(sigma_fs, kalcal_rmse[2, 0, 0, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[2].plot(sigma_fs, kalcal_rmse[2, 0, 0, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[2].plot(sigma_fs, kalcal_rmse[2, 0, 0, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[2].plot(sigma_fs, kalcal_rmse[2, 0, 0, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[2].plot(mins[2, 0, 0, 0], mins[2, 0, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[2].set_xscale("log")

# Set legend
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, facecolor='white', loc='lower center', ncol=5)

# Change x-ticks
for ax in axes:
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.4g'))

# Save plot
print("# Saving figure 1")
plt.savefig(f"plots/kalcal/rmse-all-src-1.png")

# Plot RMSE
print("# Plotting figure 2")
fig, axes = plt.subplots(2, 3, figsize=(16, 9))

# Labels
axes[0, 0].set_ylabel("RMSE")
axes[1, 0].set_ylabel("RMSE")
axes[1, 0].set_xlabel(r"Process Noise - $\sigma_f$")
axes[1, 1].set_xlabel(r"Process Noise - $\sigma_f$")
axes[1, 2].set_xlabel(r"Process Noise - $\sigma_f$")

axes[0, 0].set_title(f"KAT-7\n100% Modelled Flux")
axes[0, 1].set_title(f"VLA-B\n100% Modelled Flux")
axes[0, 2].set_title(f"MeerKAT\n100% Modelled Flux")
axes[1, 0].set_title(f"70% Modelled Flux")
axes[1, 1].set_title(f"70% Modelled Flux")
axes[1, 2].set_title(f"70% Modelled Flux")

fig.suptitle(f"Calibration over Two Sources\nRMSE between True Gains and Calibrated Gains (Kalcal over Process Noise)")

# # Plot values
sigma_fs = np.round(np.logspace(arrays["kat7"][2][0][1], arrays["kat7"][2][0][2], n_points), decimals)
axes[0, 0].plot(sigma_fs, kalcal_rmse[0, 1, 0, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[0, 0].plot(sigma_fs, kalcal_rmse[0, 1, 0, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[0, 0].plot(sigma_fs, kalcal_rmse[0, 1, 0, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[0, 0].plot(sigma_fs, kalcal_rmse[0, 1, 0, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[0, 0].plot(mins[0, 1, 0, 0], mins[0, 1, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[0, 0].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["vlab"][2][0][1], arrays["vlab"][2][0][2], n_points), decimals)
axes[0, 1].plot(sigma_fs, kalcal_rmse[1, 1, 0, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[0, 1].plot(sigma_fs, kalcal_rmse[1, 1, 0, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[0, 1].plot(sigma_fs, kalcal_rmse[1, 1, 0, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[0, 1].plot(sigma_fs, kalcal_rmse[1, 1, 0, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[0, 1].plot(mins[1, 1, 0, 0], mins[1, 1, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[0, 1].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["meerkat"][2][0][1], arrays["meerkat"][2][0][2], n_points), decimals)
axes[0, 2].plot(sigma_fs, kalcal_rmse[2, 1, 0, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[0, 2].plot(sigma_fs, kalcal_rmse[2, 1, 0, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[0, 2].plot(sigma_fs, kalcal_rmse[2, 1, 0, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[0, 2].plot(sigma_fs, kalcal_rmse[2, 1, 0, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[0, 2].plot(mins[2, 1, 0, 0], mins[2, 1, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[0, 2].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["kat7"][2][1][1], arrays["kat7"][2][1][2], n_points), decimals)
axes[1, 0].plot(sigma_fs, kalcal_rmse[0, 1, 1, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[1, 0].plot(sigma_fs, kalcal_rmse[0, 1, 1, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[1, 0].plot(sigma_fs, kalcal_rmse[0, 1, 1, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[1, 0].plot(sigma_fs, kalcal_rmse[0, 1, 1, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[1, 0].plot(mins[0, 1, 1, 0], mins[0, 1, 1, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[1, 0].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["vlab"][2][1][1], arrays["vlab"][2][1][2], n_points), decimals)
axes[1, 1].plot(sigma_fs, kalcal_rmse[1, 1, 1, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[1, 1].plot(sigma_fs, kalcal_rmse[1, 1, 1, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[1, 1].plot(sigma_fs, kalcal_rmse[1, 1, 1, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[1, 1].plot(sigma_fs, kalcal_rmse[1, 1, 1, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[1, 1].plot(mins[1, 1, 1, 0], mins[1, 1, 1, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[1, 1].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["meerkat"][2][1][1], arrays["meerkat"][2][1][2], n_points), decimals)
axes[1, 2].plot(sigma_fs, kalcal_rmse[2, 1, 1, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[1, 2].plot(sigma_fs, kalcal_rmse[2, 1, 1, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[1, 2].plot(sigma_fs, kalcal_rmse[2, 1, 1, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[1, 2].plot(sigma_fs, kalcal_rmse[2, 1, 1, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[1, 2].plot(mins[2, 1, 1, 0], mins[2, 1, 1, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[1, 2].set_xscale("log")

# Set legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, facecolor='white', loc='lower center', ncol=5)

for ax1 in axes:
    for ax in ax1:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.4g'))

# Save plot
print("# Saving figure 2")
plt.savefig(f"plots/kalcal/rmse-all-src-2.png")

# Plot RMSE
print("# Plotting figure 3")
fig, axes = plt.subplots(2, 2, figsize=(16, 9))

# Labels
axes[0, 0].set_ylabel("RMSE")
axes[1, 0].set_ylabel("RMSE")
axes[1, 0].set_xlabel(r"Process Noise - $\sigma_f$")
axes[1, 1].set_xlabel(r"Process Noise - $\sigma_f$")

axes[0, 0].set_title(f"100% Modelled Flux")
axes[0, 1].set_title(f"69% Modelled Flux")
axes[1, 0].set_title(f"47% Modelled Flux")
axes[1, 1].set_title(f"28% Modelled Flux")

fig.suptitle(f"VLA-B calibration over 100 Sources\nRMSE between True Gains and Calibrated Gains (Kalcal over Process Noise)")

# Plot values
sigma_fs = np.round(np.logspace(arrays["vlab"][100][0][1], arrays["vlab"][100][0][2], n_points), decimals)
axes[0, 0].plot(sigma_fs, kalcal_rmse[1, 2, 0, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[0, 0].plot(sigma_fs, kalcal_rmse[1, 2, 0, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[0, 0].plot(sigma_fs, kalcal_rmse[1, 2, 0, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[0, 0].plot(sigma_fs, kalcal_rmse[1, 2, 0, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[0, 0].plot(mins[1, 2, 0, 0], mins[1, 2, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[0, 0].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["vlab"][100][1][1], arrays["vlab"][100][1][2], n_points), decimals)
axes[0, 1].plot(sigma_fs, kalcal_rmse[1, 2, 1, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[0, 1].plot(sigma_fs, kalcal_rmse[1, 2, 1, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[0, 1].plot(sigma_fs, kalcal_rmse[1, 2, 1, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[0, 1].plot(sigma_fs, kalcal_rmse[1, 2, 1, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[0, 1].plot(mins[1, 2, 1, 0], mins[1, 2, 1, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[0, 1].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["vlab"][100][2][1], arrays["vlab"][100][2][2], n_points), decimals)
axes[1, 0].plot(sigma_fs, kalcal_rmse[1, 2, 2, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[1, 0].plot(sigma_fs, kalcal_rmse[1, 2, 2, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[1, 0].plot(sigma_fs, kalcal_rmse[1, 2, 2, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[1, 0].plot(sigma_fs, kalcal_rmse[1, 2, 2, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[1, 0].plot(mins[1, 2, 2, 0], mins[1, 2, 2, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[1, 0].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["vlab"][100][3][1], arrays["vlab"][100][3][2], n_points), decimals)
axes[1, 1].plot(sigma_fs, kalcal_rmse[1, 2, 3, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[1, 1].plot(sigma_fs, kalcal_rmse[1, 2, 3, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[1, 1].plot(sigma_fs, kalcal_rmse[1, 2, 3, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[1, 1].plot(sigma_fs, kalcal_rmse[1, 2, 3, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[1, 1].plot(mins[1, 2, 3, 0], mins[1, 2, 3, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[1, 1].set_xscale("log")

# Set legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, facecolor='white', loc='lower center', ncol=5)

for ax1 in axes:
    for ax in ax1:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.4g'))

# Save plot
print("# Saving figure 3")
plt.savefig(f"plots/kalcal/rmse-vlab-src-100.png")

# Plot RMSE
print("# Plotting figure 4")
fig, axes = plt.subplots(2, 2, figsize=(16, 9))

# Labels
axes[0, 0].set_ylabel("RMSE")
axes[1, 0].set_ylabel("RMSE")
axes[1, 0].set_xlabel(r"Process Noise - $\sigma_f$")
axes[1, 1].set_xlabel(r"Process Noise - $\sigma_f$")

axes[0, 0].set_title(f"100% Modelled Flux")
axes[0, 1].set_title(f"69% Modelled Flux")
axes[1, 0].set_title(f"47% Modelled Flux")
axes[1, 1].set_title(f"28% Modelled Flux")

fig.suptitle(f"MeerKAT calibration over 100 Sources\nRMSE between True Gains and Calibrated Gains (Kalcal over Process Noise)")

# Plot values
sigma_fs = np.round(np.logspace(arrays["meerkat"][100][0][1], arrays["meerkat"][100][0][2], n_points), decimals)
axes[0, 0].plot(sigma_fs, kalcal_rmse[2, 2, 0, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[0, 0].plot(sigma_fs, kalcal_rmse[2, 2, 0, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[0, 0].plot(sigma_fs, kalcal_rmse[2, 2, 0, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[0, 0].plot(sigma_fs, kalcal_rmse[2, 2, 0, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[0, 0].plot(mins[2, 2, 0, 0], mins[2, 2, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[0, 0].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["meerkat"][100][1][1], arrays["meerkat"][100][1][2], n_points), decimals)
axes[0, 1].plot(sigma_fs, kalcal_rmse[2, 2, 1, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[0, 1].plot(sigma_fs, kalcal_rmse[2, 2, 1, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[0, 1].plot(sigma_fs, kalcal_rmse[2, 2, 1, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[0, 1].plot(sigma_fs, kalcal_rmse[2, 2, 1, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[0, 1].plot(mins[2, 2, 1, 0], mins[2, 2, 1, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[0, 1].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["meerkat"][100][2][1], arrays["meerkat"][100][2][2], n_points), decimals)
axes[1, 0].plot(sigma_fs, kalcal_rmse[2, 2, 2, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[1, 0].plot(sigma_fs, kalcal_rmse[2, 2, 2, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[1, 0].plot(sigma_fs, kalcal_rmse[2, 2, 2, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[1, 0].plot(sigma_fs, kalcal_rmse[2, 2, 2, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[1, 0].plot(mins[2, 2, 2, 0], mins[2, 2, 2, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[1, 0].set_xscale("log")

sigma_fs = np.round(np.logspace(arrays["meerkat"][100][3][1], arrays["meerkat"][100][3][2], n_points), decimals)
axes[1, 1].plot(sigma_fs, kalcal_rmse[2, 2, 3, :, 0], linestyle="-", marker="x", color="tab:red", label="Filter")
axes[1, 1].plot(sigma_fs, kalcal_rmse[2, 2, 3, :, 1], linestyle="-", marker="x", color="forestgreen", label="Smoother")
# axes[1, 1].plot(sigma_fs, kalcal_rmse[2, 2, 3, :, 2], linestyle="-", marker="x", color="violet", label="Filter:RCS")
# axes[1, 1].plot(sigma_fs, kalcal_rmse[2, 2, 3, :, 3], linestyle="-", marker="x", color="limegreen", label="Smoother:RCS")
axes[1, 1].plot(mins[2, 2, 3, 0], mins[2, 2, 3, 1], "o", markersize=8, color="midnightblue", label="Minimum MSE")
axes[1, 1].set_xscale("log")

# Set legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, facecolor='white', loc='lower center', ncol=5)

for ax1 in axes:
    for ax in ax1:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.4g'))

# Save plot
print("# Saving figure 4")
plt.savefig(f"plots/kalcal/rmse-meerkat-src-100.png")

# Show plots
plt.show()
