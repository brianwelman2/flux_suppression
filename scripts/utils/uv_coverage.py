from casacore.tables import table
import matplotlib.pyplot as plt
import numpy as np

arrays = ["kat7", "vlab", "meerkat"]

for array in arrays:
    with table(f"ms/{array}.ms", ack=False) as tb:
        uvw = tb.getcol("UVW")
        _, tbin_indices, tbin_counts = np.unique(tb.getcol("TIME"), return_index=True, return_counts=True)

    n_time = len(tbin_indices)
    n_row = uvw.shape[0]
    n_bl = n_row//n_time

    baselines = n_time

    uvw = uvw[0:(n_bl * baselines)]

    plt.figure(array)
    plt.title(f"UV-Coverage for {array.upper()}")
    plt.xlabel("v")
    plt.ylabel("u")

    plt.plot(uvw[:, 0], uvw[:, 1], "+", color="red")
    plt.plot(-uvw[:, 0], -uvw[:, 1], "+", color="blue")

    xlim = np.array(plt.xlim())
    ylim = np.array(plt.ylim())

    axis_max = np.max((xlim, ylim))

    plt.xlim((-axis_max, axis_max))
    plt.ylim((-axis_max, axis_max))

plt.show()    