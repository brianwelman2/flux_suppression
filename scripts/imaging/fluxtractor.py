import os, contextlib
os.chdir("/home/welman/masters/projects/flux_suppression")

import numpy as np
import packratt
from pyrap.tables import table
from ducc0.wgridder import ms2dirty, dirty2ms
from africanus.constants import c as lightspeed
from africanus.calibration.utils import corrupt_vis, chunkify_rows
import matplotlib.pyplot as plt
from pfb.utils.fits import set_wcs
import Tigger
import astLib.astWCS

# solves Ax = b given an operator A and b
def pcg(A,
        b,
        x0,
        M=None,
        tol=1e-5,
        maxit=500,
        report_freq=10):

    if M is None:
        def M(x): return x

    r = A(x0) - b
    y = M(r)
    p = -y
    rnorm = np.vdot(r, y)
    if np.isnan(rnorm) or rnorm == 0.0:
        eps0 = 1.0
    else:
        eps0 = rnorm
    k = 0
    x = x0
    eps = 1.0
    stall_count = 0
    while eps > tol and k < maxit:
        xp = x.copy()
        rp = r.copy()
        Ap = A(p)
        rnorm = np.vdot(r, y)
        alpha = rnorm / np.vdot(p, Ap)
        x = xp + alpha * p
        r = rp + alpha * Ap
        y = M(r)
        rnorm_next = np.vdot(r, y)

        beta = rnorm_next / rnorm
        p = beta * p - y
        rnorm = rnorm_next
        k += 1
        epsx = np.linalg.norm(x - xp) / np.linalg.norm(x)
        epsn = rnorm / eps0
        epsp = eps
        eps = np.maximum(epsx, epsn)

    #     if not k % report_freq:
    #         print("At iteration %i eps = %f" % (k, eps))

    # if k >= maxit:
    #     print("Max iters reached. eps = %f." % eps)
    # else:
    #     print("Success, converged after %i iters" % k)
    return x

def calculate_restored_flux(msname, sky_model_file, gains_file, model_column,
            clean_vis_column, model_file, fluxes_file, nthreads, tol):

    ms = table(f'{msname}', ack=False)
    uvw = ms.getcol('UVW')
    time = ms.getcol('TIME')
    ant1 = ms.getcol('ANTENNA1')
    ant2 = ms.getcol('ANTENNA2')
    ms.close()
    utime = np.unique(time)

    spw = table(f'{msname}::SPECTRAL_WINDOW', ack=False)
    freq = spw.getcol('CHAN_FREQ')[0]
    nfreq = 1

    with open(model_file, "rb") as file:
        model = np.load(file)
        Ix = np.load(file)
        Iy = np.load(file)
        cell_rad = np.load(file)

    nx, ny = model.shape
    mask = np.where(model, 1.0, 0)

    with table(msname, ack=False) as tb:
        model_vis = tb.getcol(model_column)[:, None].astype(np.complex128)
        vis = tb.getcol(clean_vis_column)[:, :, 0].astype(np.complex128)


    model_vis = dirty2ms(uvw=uvw, freq=freq, dirty=model,
                         pixsize_x=cell_rad, pixsize_y=cell_rad,
                         epsilon=tol, nthreads=nthreads, do_wstacking=True)[:, None, None]

    #jones = np.load(f"gains/true/vlab/vlab-src-{sources}-gains.npy")[..., 0]
    jones = np.load(gains_file)[..., 0, 0]
    jones = jones[..., None]

    _, tbin_idx, tbin_counts = chunkify_rows(time, -1)

    # V = Jp int I kpq dl dm/n Jq.H
    # V = G R mask x   G = Mueller term,  G = Jp Jq.H,  G.H G = Jq Jq.H Jp.H Jp
    G = corrupt_vis(tbin_idx, tbin_counts, ant1, ant2,
                    jones, np.ones_like(model_vis))[:, :, 0]

    # x = (R.H G.H G R)inv R.H G.H V

    dirty = ms2dirty(uvw=uvw, freq=freq, ms=G.conj() * vis,
                     npix_x=nx, npix_y=ny,
                     pixsize_x=cell_rad, pixsize_y=cell_rad,
                     epsilon=tol, nthreads=nthreads, do_wstacking=True)

    W = (G.conj() * G).real

    def hess(x):
        tmp = dirty2ms(uvw=uvw, freq=freq, dirty=mask * x,
                       pixsize_x=cell_rad, pixsize_y=cell_rad,
                       epsilon=tol, nthreads=nthreads, do_wstacking=True)
        res = ms2dirty(uvw=uvw, freq=freq, ms=tmp, wgt=W,
                       npix_x=nx, npix_y=ny,
                       pixsize_x=cell_rad, pixsize_y=cell_rad,
                       epsilon=tol, nthreads=nthreads, do_wstacking=True)
        return mask * res

    model_rec = pcg(hess, mask * dirty, x0=np.zeros_like(dirty), tol=tol)

    restored_flux = model_rec[Ix, Iy]
    true_flux = model[Ix, Iy]

    with open(fluxes_file, "wb") as file:
        np.save(file, restored_flux)


if __name__ == "__main__":
    ms = "meerkat"
    sources = 100
    fp = 100
    sigma_f = 1e-5
    nthreads = 8
    tol = 1e-7

    msname = f"ms/{ms}.ms"
    sky_model_file = f"skymodels/{ms}/model-fp-{fp}-src-{sources}.txt"
    gains_file = f"gains/kalcal/{ms}/src-{sources}/fp-{fp}/sigma_{sigma_f}/full-smoother.npy"
    model_file = f"models/kalcal/{ms}/src-{sources}/fp-{fp}/model.npy"
    model_column = f"FP{fp}_SRC{sources}_MODEL"
    clean_vis_column = f"CLEAN_FP{fp}_SRC{sources}_DATA"
    fluxes_file = f"fluxes/kalcal/{ms}/src-{sources}/fp-{fp}/sigma_{sigma_f}/fluxes.npy"

    if not os.path.isdir(f"fluxes/kalcal/{ms}/src-{sources}/fp-{fp}/sigma_{sigma_f}/"):
        os.mkdir(f"fluxes/kalcal/{ms}/src-{sources}/fp-{fp}/sigma_{sigma_f}/")

    calculate_restored_flux(msname, sky_model_file, gains_file, model_column,
                                clean_vis_column, model_file, fluxes_file, nthreads, tol)

    with open(fluxes_file, "rb") as file:
        restored_flux = np.load(file)

    with open(model_file, "rb") as file:
        model = np.load(file)
        Ix = np.load(file)
        Iy = np.load(file)

    true_flux = model[Ix, Iy]

    # fig, ax = plt.subplots(1, 2)
    # ax[0].imshow(model_rec)
    # ax[1].imshow(model)
    # plt.show()
    # exit()
    weighted_average = np.sum(restored_flux)/np.sum(true_flux)
    abs_diff_average = np.sum(np.abs(1 - restored_flux/true_flux) * true_flux)/np.sum(true_flux)
    flux_change = restored_flux - true_flux
    # sum(abs(1 - r/t) * t)/sum(t)
    # weighted mean = sum(r)/sum(t) <-- cyndie thesis
    # flux amplification/suppression = r - t

    print(restored_flux)
    print(true_flux)
    print(weighted_average)
    print(abs_diff_average)

    from matplotlib.cm import ScalarMappable

    plt.style.use('ggplot')

    n = 100
    plt.figure(figsize=(16, 9))
    my_cmap = plt.cm.get_cmap('plasma')
    colors = my_cmap(true_flux)
    plt.bar(np.arange(1, n + 1), flux_change, label="Flux Difference", color=colors)
    plt.xlabel("Sources")
    plt.ylabel("Flux Difference")
    plt.xlim((0, n + 1))
    plt.ylim((-0.02, 0.02))
    plt.title("Flux Difference between Restored and True Flux\nMeerKAT, 100% modelled sources, sigma_f=1e-5")

    sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(np.min(true_flux), np.max(true_flux)))
    sm.set_array([])
    cbar = plt.colorbar(sm)
    cbar.set_label("Source Flux")

    # plt.plot(np.arange(1, n + 1), true_flux[:n]/true_flux[:n], "k--", label="True Flux")
    # sc = plt.scatter(np.arange(1, n + 1), restored_flux[:n]/true_flux[:n],
    #                     marker="x", s=20, c=np.log(restored_flux[:n]),
    #                     cmap=plt.get_cmap("plasma"), label="Restored Flux")

    # plt.xlabel("Sources")
    # plt.ylabel("Flux Ratio")
    # plt.ylim((0.945, 1.085))
    # plt.title("Flux Ratio of Restored vs True Flux\nVLA-B, 100% modelled sources, sigma_f=1e-5")
    # cb = plt.colorbar(sc)
    # cb.set_label("log(Flux)")

    plt.savefig("bar_true_meerkat_100.png")
    plt.show()
