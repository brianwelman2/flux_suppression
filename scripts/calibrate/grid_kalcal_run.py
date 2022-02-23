import os
os.chdir("/home/welman/masters/projects/flux_suppression")

import numpy as np
from time import time
from datetime import datetime
from kalcal.calibration.vanilla import calibrate
from tqdm import tqdm

# All time solution intervals
sigma_n = 1.0
step_control = 0.5
n_points = 16
nthreads = 24
vault_path = "/net/tina/vault2-tina/welman/flux_suppression/"
decimals = 16
start = time()

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
logfile = open(f"logs/grid-kalcal-run-{datetime.today().strftime('%d-%m-%y')}.log", "w").close()

with tqdm(total=n_points * 17) as pbar:
      for ant, sources in arrays.items():
            for source, setup in sources.items():
                  for (percent, lb, ub) in setup:
                        sigma_fs = np.round(np.logspace(lb, ub, n_points), 16)
                        for sigma_f in sigma_fs:
                              # Header
                              # print(f"!~~~~~~ TIME = {datetime.now()}, ANTENNA = {ant}, SRC = {source}, FP = {percent}, SIGMA_F = {sigma_f} ~~~~~~!")
                              with open(f"logs/grid-kalcal-run-{datetime.today().strftime('%d-%m-%y')}.log", "a") as logfile:
                                    logfile.write(f"!~~~~~~ TIME = {datetime.now()}, ANTENNA = {ant}, SRC = {source}, FP = {percent}, SIGMA_F = {sigma_f} ~~~~~~!\n")

                              # Create folder if not there
                              if not os.path.isdir(vault_path + f"gains/kalcal/{ant}/src-{source}/fp-{percent}/sigma_{sigma_f}/"):
                                    os.mkdir(vault_path + f"gains/kalcal/{ant}/src-{source}/fp-{percent}/sigma_{sigma_f}/")

                              ADJ = 0

                              if source == 1 or source == 2:
                                    MODEL_PREFIX = ""
                                    DATA_PREFIX = ""
                                    if source == 2 and int(percent) == 70:
                                          ADJ = 1
                              else:
                                    MODEL_PREFIX = f"FP{percent}_"
                                    DATA_PREFIX = f"FP100_"

                              # Kalcal calibrate
                              calibrate(
                                    f"ms/{ant}.ms",
                                    filter=1,
                                    smoother=1,
                                    algorithm="NUMBA",
                                    sigma_f=float(sigma_f),
                                    sigma_n=float(sigma_n),
                                    step_control=step_control,
                                    model_column=MODEL_PREFIX + f"SRC{source - ADJ}_MODEL",
                                    vis_column=DATA_PREFIX + f"SRC{source}_DATA",
                                    weight_column="WEIGHT",
                                    out_filter=vault_path + f"gains/kalcal/{ant}/src-{source}/fp-{percent}/sigma_{sigma_f}/full-filter.npy",
                                    out_smoother=vault_path + f"gains/kalcal/{ant}/src-{source}/fp-{percent}/sigma_{sigma_f}/full-smoother.npy",
                                    out_data="",
                                    out_weight="",
                                    ncpu=nthreads,
                                    yaml=None
                              )

                              pbar.update(1)

print(f"Done in {time() - start} s")
