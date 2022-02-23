import os
os.chdir("/home/brian/Degrees/masters/code/python/projects/flux_suppression")

import numpy as np
import subprocess as sbp
from time import time
from datetime import datetime

# All time solution intervals
tints = np.arange(12, 252, 12)
tints = np.insert(tints, 0, 8)
tints = np.insert(tints, 0, 4)
tints = np.insert(tints, 0, 1)

arrays = {
      "kat7": {
            1: [
                  (100, 100, 30, 0.2)
            ],

            2: [
                  (100, 100, 30, 0.2),
                  (70, 100, 30, 0.2)
            ]

      },
      "vlab": {
            1: [
                  (100, 100, 0.6, 0.1)
            ], 
            
            2: [
                  (100, 2000, 1, 0.1),
                  (70, 2000, 1, 0.1)
            ],

            100: [
                  (100, 2000, 2, 0.1),
                  (70, 2000, 2, 0.1),
                  (50, 2000, 2, 0.1),
                  (30, 2000, 2, 0.1)
            ]
      },
      "meerkat": {
            1: [
                  (100, 100, 0.6, 0.1)
            ],
            
            2: [
                  (100, 2000, 1, 0.1),
                  (70, 2000, 1, 0.1)
            ],

            100: [
                  (100, 2000, 2, 0.1),
                  (70, 2000, 2, 0.1),
                  (50, 2000, 2, 0.1),
                  (30, 2000, 2, 0.1)
            ]
      }
}

start = time()
logfile = open(f"logs/grid-cubical-run-{datetime.today().strftime('%d-%m-%y')}.log", "w")

for ant, sources in arrays.items():
      for source, setup in sources.items():
            for (percent, psize, pscale, mgain) in setup:
                  for tint in tints:
                        # Header
                        print(f"!~~~~~~ ANTENNA = {ant}, SRC = {source}, FP = {percent}, TINT = {tint} ~~~~~~!")
                        logfile.write(f"!~~~~~~ ANTENNA = {ant}, SRC = {source}, FP = {percent}, TINT = {tint} ~~~~~~!\n")
                        
                        # Create folder if not there
                        if not os.path.isdir(f"gains/cubical/{ant}/src-{source}/fp-{percent}/tint_{tint}/"):
                              os.mkdir(f"gains/cubical/{ant}/src-{source}/fp-{percent}/tint_{tint}/")

                        if not os.path.isdir(f"images/cubical-img/{ant}/src-{source}/fp-{percent}/tint_{tint}/"):
                              os.mkdir(f"images/cubical-img/{ant}/src-{source}/fp-{percent}/tint_{tint}/")
                        
                        ADJ = 0

                        if source == 1 or source == 2:
                              MODEL_PREFIX = ""
                              DATA_PREFIX = ""
                              if source == 2 and int(percent) == 70:
                                    ADJ = 1
                        else:
                              MODEL_PREFIX = f"FP{percent}_"
                              DATA_PREFIX = f"FP100_"

                        # Cubical Cmd Args
                        cubical_args = ["gocubical", "configs/cubical/cubical-template.parset",
                                    "--model-list", MODEL_PREFIX + f"SRC{source - ADJ}_MODEL",
                                    "--data-ms", f"ms/{ant}.ms",
                                    "--data-column", DATA_PREFIX + f"SRC{source}_DATA", 
                                    "--g-time-int", f"{tint}",
                                    "--out-dir", f"cubical/{ant}/src-{source}/fp-{percent}/tint_{tint}/",
                                    "--g-save-to", f"gains/cubical/{ant}/src-{source}/fp-{percent}/tint_{tint}/cubical.db"]

                        # Run Cubical Cmd 
                        sbp.Popen(cubical_args).wait()
                        
                        # WSClean Cmd Args
                        wsclean_args = [ "wsclean", 
                                    "-size", f"{psize}", f"{psize}", 
                                    "-scale", f"{pscale}asec",
                                    "-mgain", f"{mgain}", 
                                    "-data-column", f"CORRECTED_DATA",
                                    "-name", f"images/cubical-img/{ant}/src-{source}/fp-{percent}/tint_{tint}/{ant}-fp-{percent}-src-{source}-{tint}", 
                                    "-niter", "10000000", 
                                    "-weight", "uniform", 
                                    f"ms/{ant}.ms"]

                        # Run WSClean Cmd 
                        sbp.Popen(wsclean_args).wait() 

os.remove("cubical.last")
os.remove("montblanc.log")
logfile.close()
print(f"Done in {time() - start} s")