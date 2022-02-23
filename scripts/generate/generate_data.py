import os, contextlib
os.chdir("/home/brian/Degrees/masters/code/python/projects/flux_suppression")

import numpy as np
from pyrap.tables import table, maketabdesc, makearrcoldesc
from ducc0.wgridder import ms2dirty, dirty2ms
from africanus.constants import c as lightspeed
from africanus.calibration.utils import corrupt_vis, chunkify_rows
from astropy.io import fits
from astropy import wcs
import astLib.astWCS
import Tigger
import subprocess as sbp
from pfb.utils.fits import set_wcs
from tqdm import tqdm


def gen_data(ms_name, sky_model, gains_file, sigma_n, model_column, data_column):
    with table(ms_name, ack=False) as tb:
        uvw = tb.getcol('UVW')
        time = tb.getcol('TIME')
        ant1 = tb.getcol('ANTENNA1')
        ant2 = tb.getcol('ANTENNA2')

    n_row = uvw.shape[0]
    u_time = np.unique(time)
    n_time = u_time.size
    n_ant = int(np.maximum(ant1.max(), ant2.max()) + 1)

    with table(ms_name + "::SPECTRAL_WINDOW", ack=False) as tb:
        freq = tb.getcol('CHAN_FREQ')[0]

    n_chan = freq.size

    with table(ms_name + "::FIELD", ack=False) as tb:
        radec = tb.getcol('PHASE_DIR')[0][0]

    u_max = np.abs(uvw[:, 0]).max()
    v_max = np.abs(uvw[:, 1]).max()
    uv_max = np.maximum(u_max, v_max)

    # twice Nyquist freq
    cell = 0.5/(2*uv_max*freq[0]/lightspeed)

    nx = ny = 512

    header = set_wcs(cell*3600, cell*3600, nx, ny, radec, freq)

    smodel = Tigger.load(sky_model, verbose=False)    
    wcs1 = astLib.astWCS.WCS(header, mode="pyfits")
    model = np.zeros((nx, ny), dtype=np.float64)
    Ix, Iy = [], []
    for src in smodel.sources:
        y, x = wcs1.wcs2pix(np.rad2deg(src.pos.ra), np.rad2deg(src.pos.dec))
        x, y = int(x), int(y)
        Ix.append(x)
        Iy.append(y)
        
        model[x, y] = src.flux.I

    Ix = np.array(Ix)
    Iy = np.array(Iy)

    model_vis = dirty2ms(uvw=uvw, freq=freq, dirty=model,
                         pixsize_x=cell, pixsize_y=cell,
                         epsilon=1e-7, nthreads=8, do_wstacking=True)[:, None, None]

    jones = np.load(gains_file)[..., 0]
    jones = jones[..., None]

    _, tbin_idx, tbin_counts = chunkify_rows(time, -1)

     
    clean_vis = corrupt_vis(tbin_idx, tbin_counts, ant1, ant2,
                                jones, model_vis)

    noise = np.random.normal(loc=0.0, scale=sigma_n, size=clean_vis.shape) \
            + 1.0j * np.random.normal(loc=0.0, scale=sigma_n, size=clean_vis.shape)

    vis = clean_vis + noise

    npad = ((0, 0), (0, 0), (0, 3))
    clean_vis = np.pad(clean_vis, pad_width=npad, mode='constant', constant_values=0j).astype(np.complex128)
    vis = np.pad(vis, pad_width=npad, mode='constant', constant_values=0j).astype(np.complex128)
    model_vis = np.pad(model_vis[..., 0], pad_width=npad, mode='constant', constant_values=0j).astype(np.complex128)

    with table(ms_name, ack=False, readonly=False) as tb:
        datadmi = tb.getdminfo('DATA')   
        datadmi["NAME"] = data_column
        tb.addcols(maketabdesc(makearrcoldesc(datadmi["NAME"], 0.0j, ndim=2)), datadmi)
        datadmi["NAME"] = "CLEAN_" + data_column
        tb.addcols(maketabdesc(makearrcoldesc(datadmi["NAME"], 0.0j, ndim=2)), datadmi)
        datadmi["NAME"] = model_column
        tb.addcols(maketabdesc(makearrcoldesc(datadmi["NAME"], 0.0j, ndim=2)), datadmi)

        tb.putcol(data_column, vis)
        tb.putcol("CLEAN_" + data_column, clean_vis)
        tb.putcol(model_column, model_vis)

def main():

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
        },
    }

    sigma_n = np.sqrt(2)

    for ant, sources in arrays.items():
            for source, percents in sources.items():
                for percent in percents:            
                    print(f":: ANT={ant}\tSRC={source}\tPRT={percent}")                    
                    if source == 1:
                        model_column = "SRC1_MODEL"
                        data_column = "SRC1_DATA"
                        sky_model = f"skymodels/{ant}/model-src-1.txt"
                    elif source == 2 and percent == 100:
                        model_column = "SRC2_MODEL"
                        data_column = "SRC2_DATA"
                        sky_model = f"skymodels/{ant}/model-src-2.txt"
                    elif source == 2 and percent == 70:
                        continue
                    else:
                        if ant == "kat7":
                            continue
                        model_column = f"FP{percent}_SRC100_MODEL"
                        data_column = f"FP{percent}_SRC100_DATA"
                        sky_model = f"skymodels/{ant}/model-fp-{percent}-src-100.txt"
                            
                    # Calculated restored flux
                    ms_name = f"ms/{ant}.ms"
                    gains_file = f"gains/true/{ant}/{ant}-src-{source}-gains.npy"
                    
                    gen_data(ms_name, sky_model, gains_file, sigma_n, model_column, data_column)

if __name__ == "__main__":
    main()
    print("Done")    