import os
os.chdir("/home/welman/masters/projects/flux_suppression")

import numpy as np
from time import time
from datetime import datetime
from kalcal.calibration.vanilla import calibrate
from tqdm import tqdm
from fluxtractor import calculate_restored_flux
from astropy.io import fits

# All time solution intervals
sigma_n = 1.0
step_control = 0.5
n_points = 16
tol = 1e-7
vault_path = "/net/tina/vault2-tina/welman/flux_suppression/"
nthreads = 24
start = time()

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

logfile = open(f"logs/grid-image-run-{datetime.today().strftime('%d-%m-%y')}.log", "w").close()

with tqdm(total=n_points * 17) as pbar:
    for ant, sources in arrays.items():
        for source, setup in sources.items():
            for (percent, lb, ub) in setup:
                sigma_fs = np.round(np.logspace(lb, ub, n_points), 7)
                for sigma_f in sigma_fs:
                # Header
                # print(f"!~~~~~~ TIME = {datetime.now()}, ANTENNA = {ant}, SRC = {source}, FP = {percent}, SIGMA_F = {sigma_f} ~~~~~~!")
                    with open(f"logs/grid-image-run-{datetime.today().strftime('%d-%m-%y')}.log", "a") as logfile:
                        logfile.write(f"!~~~~~~ TIME = {datetime.now()}, ANTENNA = {ant}, SRC = {source}, FP = {percent}, SIGMA_F = {sigma_f} ~~~~~~!\n")

                     # Create folder if not there
                    if not os.path.isdir(vault_path + f"fluxes/kalcal/{ant}/src-{source}/fp-{percent}/sigma_{sigma_f}/"):
                        os.mkdir(vault_path + f"fluxes/kalcal/{ant}/src-{source}/fp-{percent}/sigma_{sigma_f}/")

                    if source == 1 or (source == 2 and percent == 70):
                        data_column = "SRC1_DATA"
                        model_column = "SRC1_MODEL"
                        sky_model = f"skymodels/{ant}/model-src-1.txt"
                    elif source == 2 and percent == 100:
                        data_column = "SRC2_DATA"
                        model_column = "SRC2_MODEL"
                        sky_model = f"skymodels/{ant}/model-src-2.txt"
                    else:
                        data_column = f"FP{percent}_SRC100_DATA"
                        model_column = f"FP{percent}_SRC100_MODEL"
                        sky_model = f"skymodels/{ant}/model-fp-{percent}-src-100.txt"

                    # Calculated restored flux
                    msname = f"ms/{ant}.ms"
                    sky_model_file = f"skymodels/{ant}/model-fp-{percent}-src-{source}.txt"
                    gains_file = vault_path + f"gains/kalcal/{ant}/src-{source}/fp-{percent}/sigma_{sigma_f}/full-smoother.npy"
                    model_file = f"models/kalcal/{ant}/src-{source}/fp-{percent}/model.npy"
                    fluxes_file = vault_path + f"fluxes/kalcal/{ant}/src-{source}/fp-{percent}/sigma_{sigma_f}/fluxes.npy"

                    calculate_restored_flux(msname, sky_model_file, gains_file,
                                            model_column, "CLEAN_" + data_column,
                                            model_file, fluxes_file, nthreads, tol)

                    pbar.update(1)

print(f"Done in {time() - start} s")
