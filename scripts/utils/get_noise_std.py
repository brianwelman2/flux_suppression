import os
os.chdir("/home/brian/Degrees/masters/code/python/projects/flux_suppression")

from casacore.tables import table
import numpy as np

msname = "ms/kalcal.ms"

with table(msname) as tb:
    noise = tb.getcol("NOISE")[:, :, (0, 3)]

gains = np.load("gains/true_gains.npy")[..., (0, 3)]
gdiff = (gains[1:] - gains[0:-1]).flatten()

var_n = np.var(noise)
sigma_n = np.sqrt(var_n/2)

var_f = np.var(gdiff)
sigma_f = np.sqrt(var_f/2)

print(f"MEASURE: complex_var = {var_n}, axis_std (sigma_n) = {sigma_n}")
print(f"STATE: complex_var = {var_f}, axis_std (sigma_f) = {sigma_f}")