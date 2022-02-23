import os, contextlib
os.chdir("/home/welman/masters/projects/flux_suppression")

import numpy as np
from ducc0.wgridder import ms2dirty, dirty2ms
from pyrap.tables import table
import Tigger
from africanus.constants import c as lightspeed
from africanus.calibration.utils import corrupt_vis, chunkify_rows
from africanus.coordinates import radec_to_lm, lm_to_radec
from ducc0.fft import good_size
from tqdm import tqdm
from datetime import datetime

# Params
ms = "vlab"
sources = 100
fp = 100
sigma_f = 1e-5
nthreads = 8
tol = 1e-7
sigma_n = 1.0

# Paths
msname = f"ms/{ms}.ms"
sky_model_file = f"skymodels/{ms}/model-fp-{fp}-src-{sources}.txt"
gains_file = f"gains/true/{ms}/{ms}-src-{sources}-gains.npy"
model_column = f"FP{fp}_SRC{sources}_MODEL"
vis_column = f"FP{fp}_SRC{sources}_DATA"
clean_vis_column = f"CLEAN_FP{fp}_SRC{sources}_DATA"
fluxes_file = f"fluxes/kalcal/{ms}/src-{sources}/fp-{fp}/sigma_{sigma_f}/fluxes.npy"

def simulate_data(msname, sky_model_file, gains_file, model_file,
                    model_column, vis_column, sigma_n=1.0, fov=1.0):

    # Main table
    ms = table(f'{msname}', ack=False)
    uvw = ms.getcol('UVW')
    time = ms.getcol('TIME')
    ant1 = ms.getcol('ANTENNA1')
    ant2 = ms.getcol('ANTENNA2')
    ms.close()
    utime = np.unique(time)

    # Spectral table
    spw = table(f'{msname}::SPECTRAL_WINDOW', ack=False)
    freq = spw.getcol('CHAN_FREQ')[0]
    nfreq = 1

    # UVmax
    umax = np.abs(uvw[:, 0]).max()
    vmax = np.abs(uvw[:, 1]).max()
    uvmax = np.maximum(umax, vmax)

    # twice Nyquist freq
    cell_rad = 0.9/(2*uvmax*freq[0]/lightspeed)
    cell_deg = np.rad2deg(cell_rad)
    cell_asec = cell_deg * 3600


    # Select good pixel size for fft
    npix = good_size(int(np.ceil(fov/cell_deg)))

    while npix % 2:
        npix += 1
        npix = good_size(npix)

    # Set pixel size
    nx = npix
    ny = npix

    # Models
    model = np.zeros((nx, ny), dtype=np.float64)
    smodel = Tigger.load(sky_model_file, verbose=False)

    # LM Grid
    x = -(nx//2) + np.arange(nx)
    y = -(ny//2) + np.arange(nx)
    l = np.require(x, np.float64)
    m = np.require(y, np.float64)

    # set cell_rad in radians
    l *= cell_rad
    m *= cell_rad

    # Phase Centre
    field = table(f'{msname}::FIELD', ack=False)
    phase_centre = field.getcol('PHASE_DIR').squeeze()

    # Source Coordinates
    Ix, Iy, = [], []

    # RADEC to Grid points
    for src in smodel.sources:
        radec = np.array([[src.pos.ra, src.pos.dec]])
        l0, m0 = radec_to_lm(radec, phase_centre)[0]
        dell = l - l0
        delm = m - m0

        xi = np.abs(l - l0).argmin()
        yi = np.abs(m - m0).argmin()
        model[xi, yi] = src.flux.I

        Ix.append(xi)
        Iy.append(yi)

    with open(model_file, "wb") as file:
        np.save(file, model)
        np.save(file, Ix)
        np.save(file, Iy)
        np.save(file, cell_rad)

    # Create model visibilities
    model_vis = dirty2ms(uvw=uvw, freq=freq, dirty=model,
                         pixsize_x=cell_rad, pixsize_y=cell_rad,
                         epsilon=1e-7, nthreads=8, do_wstacking=True)[:, None, None]

    # Load true gains
    jones = np.load(gains_file)[..., 0]
    jones = jones[..., None]

    # Time bin params
    _, tbin_idx, tbin_counts = chunkify_rows(time, -1)

    # Corrupt model visibilities with true gains
    clean_vis = corrupt_vis(tbin_idx, tbin_counts, ant1, ant2,
                                jones, model_vis)

    # Measure noise for visibilities
    noise = np.random.normal(loc=0.0, scale=sigma_n, size=clean_vis.shape)\
            + 1.0j * np.random.normal(loc=0.0, scale=sigma_n, size=clean_vis.shape)

    # Visibilities
    vis = clean_vis + noise

    # Write visibilities to table
    with table(msname, ack=False, readonly=False) as tb:
        model_desc = tb.getcoldesc("MODEL_DATA")
        model_desc["name"] = model_column
        try:
            tb.addcols(model_desc)
        except:
            pass

        vis_desc = tb.getcoldesc("DATA")
        vis_desc["name"] = vis_column
        try:
            tb.addcols(vis_desc)
        except:
            pass

        vis_desc["name"] = "CLEAN_" + vis_column
        try:
            tb.addcols(vis_desc)
        except:
            pass

        tb.putcol(model_column, model_vis[..., 0])
        tb.putcol("CLEAN_" + vis_column, clean_vis)
        tb.putcol(vis_column, clean_vis)

if __name__ == "__main__":
    arrays = {
      "kat7": {
            1: [
                  (100, 0.5)
            ],

            2: [
                  (100, 0.5),
                  (70, 0.5)
            ]

      },
      "vlab": {
            1: [
                  (100, 0.5)
            ],

            2: [
                  (100, 0.5),
                  (70, 0.5)
            ],

            100: [
                  (100, 1.0),
                  (70, 1.0),
                  (50, 1.0),
                  (30, 1.0)
            ]
      },
      "meerkat": {
            1: [
                  (100, 0.5)
            ],

            2: [
                  (100, 0.5),
                  (70, 0.5)
            ],

            100: [
                  (100, 1.0),
                  (70, 1.0),
                  (50, 1.0),
                  (30, 1.0)
            ]
      }
}

with tqdm(total=17) as pbar:
    for ant, sources in arrays.items():
        for source, setup in sources.items():
                for (percent, fov) in setup:
                    # Header
                    with open(f"logs/generate-data-{datetime.today().strftime('%d-%m-%y')}.log", "a") as logfile:
                        logfile.write(f"!~~~~~~ TIME = {datetime.now()}, ANTENNA = {ant}, SRC = {source}, FP = {percent} ~~~~~~!\n")

                    # Create folder if not there
                    if not os.path.isdir(f"models/kalcal/{ant}/src-{source}/fp-{percent}/"):
                        os.mkdir(f"models/kalcal/{ant}/src-{source}/fp-{percent}/")

                    if source == 1 or (source == 2 and percent == 70):
                        model_column = "SRC1_MODEL"
                        vis_column = "SRC1_DATA"
                        sky_model = f"skymodels/{ant}/model-src-1.txt"
                    elif source == 2 and percent == 100:
                        model_column = "SRC2_MODEL"
                        vis_column = "SRC2_DATA"
                        sky_model = f"skymodels/{ant}/model-src-2.txt"
                    else:
                        model_column = f"FP{percent}_SRC100_MODEL"
                        vis_column = f"FP{percent}_SRC100_DATA"
                        sky_model = f"skymodels/{ant}/model-fp-{percent}-src-100.txt"

                    msname = f"ms/{ant}.ms"
                    gains_file = f"gains/true/{ant}/{ant}-src-{source}-gains.npy"
                    model_file = f"models/kalcal/{ant}/src-{source}/fp-{percent}/model.npy"

                    simulate_data(msname, sky_model_file, gains_file, model_file,
                                    model_column, vis_column, sigma_n=1.0, fov=1.0)

                    pbar.update(1)
