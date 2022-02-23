import os
os.chdir("/home/brian/Degrees/masters/code/python/projects/flux_suppression")

import numpy as np
import subprocess as sbp


# Metric Function - Root Mean Square Error
def rmse(A, B):
    res = A - B
    return np.sqrt(np.vdot(res, res).real/res.size)

# Number of smoother runs to try to
n_smoothers = 6

true_gains = np.load("gains/true_gains.npy")[..., (0, 3)]

kalcal_rmse = np.zeros((2, n_smoothers), dtype=np.float64)

for i, s in enumerate(range(1, n_smoothers + 1)):
    full_gains = np.load(f"gains/kalcal/other/full-0.003-{s}-smoother.npy")[..., (0, 3), 0]
    half_gains = np.load(f"gains/kalcal/other/half-0.003-{s}-smoother.npy")[..., (0, 3), 0]

    kalcal_rmse[0, i] = rmse(true_gains, full_gains)
    kalcal_rmse[1, i] = rmse(true_gains, half_gains)

# Print for now
for i, val in enumerate(kalcal_rmse[0]): 
    print(f"Kalcal (FULL) with smoother={i + 1}: RMSE={val}")

print()

for i, val in enumerate(kalcal_rmse[1]): 
    print(f"Kalcal (HALF) with smoother={i + 1}: RMSE={val}")