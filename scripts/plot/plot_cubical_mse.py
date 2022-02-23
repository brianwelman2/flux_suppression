from casacore.tables import table
import os
os.chdir("/home/brian/Degrees/masters/code/python/projects/flux_suppression")

import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as plt
from tqdm import tqdm

# Params for plots
rcParams.update({'font.size': 12, 'font.family': 'serif'})
plt.style.use('ggplot')

# Metric Function - Root Mean Square Error
def rmse(A, B):
    res = A - B
    return np.sqrt(np.vdot(res, res).real/res.size)


arrays = {
      "kat7": {
            1: [
                  (100, 0.0002, 0.02)
            ],

            2: [
                  (100, 0.0002, 0.02),
                  (70, 0.0002, 0.02)
            ]

      },
      "vlab": {
            1: [
                  (100, 0.0002, 0.02,)
            ], 
            
            2: [
                  (100, 0.0002, 0.02),
                  (70, 0.0002, 0.02)
            ],

            100: [
                  (100, 0.0002, 0.02),
                  (70, 0.0002, 0.02),
                  (50, 0.0002, 0.02),
                  (30, 0.0002, 0.02)
            ]
      },
      "meerkat": {
            1: [
                  (100, 0.0002, 0.02,)
            ],
            
            2: [
                  (100, 0.0002, 0.02),
                  (70, 0.0002, 0.02)
            ],

            100: [
                  (100, 0.0002, 0.02),
                  (70, 0.0002, 0.02),
                  (50, 0.0002, 0.02),
                  (30, 0.0002, 0.02)
            ]
      }
}

# All time solution intervals
tints = np.arange(12, 252, 12)
tints = np.insert(tints, 0, 8)
tints = np.insert(tints, 0, 4)
tints = np.insert(tints, 0, 1)

# Minimums
mins = np.zeros((3, 3, 4, 2), dtype=np.float64)
cubical_rmse = np.zeros((3, 3, 4, len(tints)), dtype=np.float64)

with tqdm(total=391) as pbar:
      for i, (ant, sources) in enumerate(arrays.items()):
            for j, (source, setup) in enumerate(sources.items()):
                  for k, (percent, lb, ub) in enumerate(setup):
                        true_gains = np.load(f"gains/true/{ant}/{ant}-src-{source}-gains.npy")[..., (0, 3)]

                        for l, tint in enumerate(tints):
                              cubical_path = f"gains/cubical/{ant}/src-{source}/fp-{percent}/tint_{tint}/cubical.npy"
                              cubical_gains = np.load(cubical_path)[..., (0, 3)]
                              cubical_rmse[i, j, k, l] = rmse(true_gains, cubical_gains)  

                              pbar.update(1)              
                        
                        mins[i, j, k, 0] = tints[cubical_rmse[i, j, k].argmin()]
                        mins[i, j, k, 1] = cubical_rmse[i, j, k].min()

# Plot RMSE
fig, axes = plt.subplots(1, 3, figsize=(16, 9))

# Labels
axes[0].set_ylabel("RMSE")
axes[0].set_xlabel(r"Time Intervals - $\delta t$")
axes[1].set_xlabel(r"Time Intervals - $\delta t$")
axes[2].set_xlabel(r"Time Intervals - $\delta t$")

axes[0].set_title(f"KAT-7")
axes[1].set_title(f"VLA-B")
axes[2].set_title(f"MeerKAT")

fig.suptitle(f"Calibration over Single Source\nRMSE between True Gains and Calibrated Gains (Cubical over Time-Intervals)")

# Plot values
axes[0].plot(tints, cubical_rmse[0, 0, 0], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[0].plot(mins[0, 0, 0, 0], mins[0, 0, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

axes[1].plot(tints, cubical_rmse[1, 0, 0], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[1].plot(mins[1, 0, 0, 0], mins[1, 0, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

axes[2].plot(tints, cubical_rmse[2, 0, 0], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[2].plot(mins[2, 0, 0, 0], mins[2, 0, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

# Set legend
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, facecolor='white', loc='lower center', ncol=3)

# Plot RMSE
fig, axes = plt.subplots(2, 3, figsize=(16, 9))

# Labels
axes[0, 0].set_ylabel("RMSE")
axes[1, 0].set_ylabel("RMSE")
axes[1, 0].set_xlabel(r"Time Intervals - $\delta t$")
axes[1, 1].set_xlabel(r"Time Intervals - $\delta t$")
axes[1, 2].set_xlabel(r"Time Intervals - $\delta t$")

axes[0, 0].set_title(f"KAT-7\n100% Modelled Flux")
axes[0, 1].set_title(f"VLA-B\n100% Modelled Flux")
axes[0, 2].set_title(f"MeerKAT\n100% Modelled Flux")
axes[1, 0].set_title(f"70% Modelled Flux")
axes[1, 1].set_title(f"70% Modelled Flux")
axes[1, 2].set_title(f"70% Modelled Flux")

fig.suptitle(f"Calibration over Two Sources\nRMSE between True Gains and Calibrated Gains (Cubical over Time-Intervals)")

# Plot values
axes[0, 0].plot(tints, cubical_rmse[0, 1, 0], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[0, 0].plot(mins[0, 1, 0, 0], mins[0, 1, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

axes[0, 1].plot(tints, cubical_rmse[1, 1, 0], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[0, 1].plot(mins[1, 1, 0, 0], mins[1, 1, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

axes[0, 2].plot(tints, cubical_rmse[2, 1, 0], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[0, 2].plot(mins[2, 1, 0, 0], mins[2, 1, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

axes[1, 0].plot(tints, cubical_rmse[0, 1, 1], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[1, 0].plot(mins[0, 1, 1, 0], mins[0, 1, 1, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

axes[1, 1].plot(tints, cubical_rmse[1, 1, 1], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[1, 1].plot(mins[1, 1, 1, 0], mins[1, 1, 1, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

axes[1, 2].plot(tints, cubical_rmse[2, 1, 1], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[1, 2].plot(mins[2, 1, 1, 0], mins[2, 1, 1, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

# Set legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, facecolor='white', loc='lower center', ncol=3)

# Plot RMSE
fig, axes = plt.subplots(2, 2, figsize=(16, 9))

# Labels
axes[0, 0].set_ylabel("RMSE")
axes[1, 0].set_ylabel("RMSE")
axes[1, 0].set_xlabel(r"Time Intervals - $\delta t$")
axes[1, 1].set_xlabel(r"Time Intervals - $\delta t$")

axes[0, 0].set_title(f"100% Modelled Flux")
axes[0, 1].set_title(f"69% Modelled Flux")
axes[1, 0].set_title(f"47% Modelled Flux")
axes[1, 1].set_title(f"28% Modelled Flux")

fig.suptitle(f"VLA-B calibration over 100 Sources\nRMSE between True Gains and Calibrated Gains (Cubical over Time-Intervals)")

# Plot values
axes[0, 0].plot(tints, cubical_rmse[1, 2, 0], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[0, 0].plot(mins[1, 2, 0, 0], mins[1, 2, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

axes[0, 1].plot(tints, cubical_rmse[1, 2, 1], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[0, 1].plot(mins[1, 2, 1, 0], mins[1, 2, 1, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

axes[1, 0].plot(tints, cubical_rmse[1, 2, 2], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[1, 0].plot(mins[1, 2, 2, 0], mins[1, 2, 2, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

axes[1, 1].plot(tints, cubical_rmse[1, 2, 3], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[1, 1].plot(mins[1, 2, 3, 0], mins[1, 2, 3, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

# Set legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, facecolor='white', loc='lower center', ncol=2)

# Plot RMSE
fig, axes = plt.subplots(2, 2, figsize=(16, 9))

# Labels
axes[0, 0].set_ylabel("RMSE")
axes[1, 0].set_ylabel("RMSE")
axes[1, 0].set_xlabel(r"Time Intervals - $\delta t$")
axes[1, 1].set_xlabel(r"Time Intervals - $\delta t$")

axes[0, 0].set_title(f"100% Modelled Flux")
axes[0, 1].set_title(f"69% Modelled Flux")
axes[1, 0].set_title(f"47% Modelled Flux")
axes[1, 1].set_title(f"28% Modelled Flux")

fig.suptitle(f"MeerKAT calibration over 100 Sources\nRMSE between True Gains and Calibrated Gains (Cubical over Time-Intervals)")

# Plot values
axes[0, 0].plot(tints, cubical_rmse[2, 2, 0], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[0, 0].plot(mins[2, 2, 0, 0], mins[2, 2, 0, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

axes[0, 1].plot(tints, cubical_rmse[2, 2, 1], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[0, 1].plot(mins[2, 2, 1, 0], mins[2, 2, 1, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

axes[1, 0].plot(tints, cubical_rmse[2, 2, 2], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[1, 0].plot(mins[2, 2, 2, 0], mins[2, 2, 2, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

axes[1, 1].plot(tints, cubical_rmse[2, 2, 3], linestyle="-", marker="x", color="tab:blue", label="Cubical")
axes[1, 1].plot(mins[2, 2, 3, 0], mins[2, 2, 3, 1], "o", markersize=8, color="midnightblue", label="Minimum RMSE")

# Set legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, facecolor='white', loc='lower center', ncol=2)

# Save
#plt.savefig(f"plots/kalcal/rmse-kalcal-over-sigma_f-{ant}-{src}.png")

# Show plots 
plt.show()