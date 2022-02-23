import os, contextlib

from tqdm.std import tqdm
os.chdir("/home/brian/Degrees/masters/code/python/projects/flux_suppression")

import cubical.param_db as pdb
from cubical.tools import logger
log = logger.getLogger("param_db")
log.verbosity(-1)

import numpy as np

# All time solution intervals
tints = np.arange(12, 252, 12)
tints = np.insert(tints, 0, 8)
tints = np.insert(tints, 0, 4)
tints = np.insert(tints, 0, 1)

arrays = {
      "kat7": {
            1: [100],            
            2: [100, 70]
      },
      "vlab": {
            1: [100],            
            2: [100, 70],
            100: [100, 70, 50, 30]
      },
      "meerkat": {
            1: [100],            
            2: [100, 70],
            100: [100, 70, 50, 30]
      }
}


for ant, sources in tqdm(arrays.items()):
    for source, percents in sources.items():
        true_gains = np.load(f"gains/true/{ant}/{ant}-src-{source}-gains.npy")
        n_time, n_ant, n_chan, n_dir, n_corr = true_gains.shape
        fint = 1

        for percent in percents:
            for tint in tints:
                with open(os.devnull, 'w') as devnull:
                    with contextlib.redirect_stderr(devnull):
                        cubical_gains = pdb.load(f"gains/cubical/{ant}/src-{source}/fp-{percent}/tint_{tint}/cubical.db")["G:gain"].get_cube()
                codex_gains = np.zeros_like(true_gains)
                for d in range(n_dir):
                    for nu in range(n_chan):
                        rc = nu // fint
                        for t in range(n_time):
                            rr = t // tint
                            for a in range(n_ant):
                                for i, c in enumerate([0, 3]):
                                    codex_gains[t, a, nu, d, c] = cubical_gains[d, rr , rc , a, i, i]

                with open(f"gains/cubical/{ant}/src-{source}/fp-{percent}/tint_{tint}/cubical.npy", "wb") as file:
                    np.save(file, codex_gains)