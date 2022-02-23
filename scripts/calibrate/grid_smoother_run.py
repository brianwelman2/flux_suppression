import os
os.chdir("/home/brian/Degrees/masters/code/python/projects/flux_suppression")

import numpy as np
import subprocess as sbp

# Number of smoother runs to try to
n_smoothers = 6

for i in range(1, n_smoothers + 1):
      # Header
      print(f"!~~~~~~ SMOOTHER = {i} ~~~~~~!")

      # Kalcal Cmd Args
      kalcal_args = [
            "kal-calibrate", "vanilla",
            "--filter", "1",
            "--smoother", f"{i}",
            "--algorithm", "NUMBA",
            "--sigma-f", f"0.003",
            "--sigma-n", f"1.0",
            "--step-control", f"0.5",
            "--model-column", "FULL_MODEL",
            "--vis-column", "FULL_DATA",
            "--weight-column", "WEIGHT",
            "--out-filter", f"gains/kalcal/other/full-0.003-{i}-filter.npy",
            "--out-smoother", f"gains/kalcal/other/full-0.003-{i}-smoother.npy",
            "--out-data", "",
            "ms/kalcal.ms"
      ]

      # Run Kalcal Cmd 
      sbp.Popen(kalcal_args).wait()

      # Kalcal Cmd Args
      kalcal_args = [
            "kal-calibrate", "vanilla",
            "--filter", "1",
            "--smoother", f"{i}",
            "--algorithm", "NUMBA",
            "--sigma-f", f"0.003",
            "--sigma-n", f"1.0",
            "--step-control", f"0.5",
            "--model-column", "HALF_MODEL",
            "--vis-column", "FULL_DATA",
            "--weight-column", "WEIGHT",
            "--out-filter", f"gains/kalcal/other/half-0.003-{i}-filter.npy",
            "--out-smoother", f"gains/kalcal/other/half-0.003-{i}-smoother.npy",
            "--out-data", "",
            "ms/kalcal.ms"
      ]

      # Run Kalcal Cmd 
      sbp.Popen(kalcal_args).wait()