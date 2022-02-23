import os
os.chdir("/home/brian/Degrees/masters/code/python/projects/flux_suppression")

from casacore.tables import table
from scipy.constants import c
import numpy as np

msname = ["kat7", "meerkat", "vlab"]

for ms in msname:
    with table(f"ms/{ms}.ms") as tb:
        uv = tb.getcol("UVW")[:, (0, 1)]

    with table(f"ms/{ms}.ms" + "::SPECTRAL_WINDOW") as tb:
        nu = tb.getcol("CHAN_FREQ").flatten()[0]

    lam = c/nu
    uv_max = np.max(np.abs(uv*lam))
    cell_size = 1/(2*uv_max)

    print(f"{ms} = {(cell_size) * 3600} asec")