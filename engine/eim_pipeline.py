"""
eim_pipeline.py  --  shared matched-simulation pipeline for the EIM/CMB parity engine.

HARDENING RATIONALE (see PREDECLARATION_AND_RULINGS.md sec.2):
    Every predeclared statistic must be compared to a MATCHED LCDM null carried through the
    *identical* processing chain -- beam, mask, homogeneous noise, and inpainting -- never to
    ideal full-sky isotropy. Mask coupling in particular biases low-l parity and axis
    statistics; a full-sky chi^2 null does not see it. This module makes the null and the
    data share one pipeline so that bias is captured by the null distribution itself.

    Nothing in this module is EIM-specific. It is standard low-l CMB machinery. The EIM content
    lives entirely in the CHOICE of statistics {g_L, alpha, H5} and the FROZEN-AXIS protocol,
    both fixed in the predeclaration.

Design notes / honest limits:
  * Inpainting is the simple iterative harmonic ("cooling" / Gerchberg-Saxton) scheme: it is
    applied IDENTICALLY to null and data, so its residual bias is absorbed by the null. It is
    NOT a constrained Gaussian realization; for a real-data submission, swap in the map-maker's
    own inpainting (e.g. the Planck common-mask diffusive inpainter) -- the API here is a
    drop-in (observe(...) returns recovered alm).
  * The bundled mask is a synthetic Galactic-cut+apodization stand-in. To run on real data,
    replace make_mask() with a load of the actual analysis mask; the rest is unchanged.
"""
import os
import numpy as np
import healpy as hp

# ------------------------------------------------------------------ config
NSIDE = 64
LMAX  = 32
NPIX  = hp.nside2npix(NSIDE)
VEC   = np.array(hp.pix2vec(NSIDE, np.arange(NPIX)))     # 3 x npix unit vectors
_CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_cl_cache.npy")


# ------------------------------------------------------------------ theory Cl
def get_cl(lmax=LMAX):
    """LCDM TT raw C_l (muK^2), C_l[0:2]=0. Cached to disk. CAMB with SW-plateau fallback."""
    if os.path.exists(_CACHE):
        cl = np.load(_CACHE)
        if len(cl) >= lmax + 1:
            cl = cl[:lmax + 1].copy(); cl[:2] = 0.0
            return cl
    try:
        import camb
        pars = camb.set_params(H0=67.36, ombh2=0.02237, omch2=0.1200,
                               ns=0.9649, As=2.1e-9, tau=0.0544, mnu=0.06)
        pars.set_for_lmax(lmax + 60)
        res = camb.get_results(pars)
        cl = np.asarray(res.get_cmb_power_spectra(pars, CMB_unit='muK',
                                                  raw_cl=True)['total'][:lmax + 1, 0])
    except Exception:
        l = np.arange(lmax + 1); cl = np.zeros(lmax + 1)
        cl[2:] = 2 * np.pi * 1000.0 / (l[2:] * (l[2:] + 1))
    cl = np.asarray(cl, float); cl[:2] = 0.0
    try:
        np.save(_CACHE, cl)
    except Exception:
        pass
    return cl


# ------------------------------------------------------------------ alm helpers
def lset_map(alm, lset, lmax=LMAX):
    a = np.zeros_like(alm)
    for L in lset:
        for m in range(L + 1):
            i = hp.Alm.getidx(lmax, L, m); a[i] = alm[i]
    return hp.alm2map(a, NSIDE)

def band_map(alm, lo, hi, lmax=LMAX):
    a = np.zeros_like(alm)
    for L in range(lo, hi + 1):
        for m in range(L + 1):
            i = hp.Alm.getidx(lmax, L, m); a[i] = alm[i]
    return hp.alm2map(a, NSIDE)


# ------------------------------------------------------------------ axis estimator
def axis_of(alm, lset, lmax=LMAX):
    """Convention-free preferred (planar) axis via power tensor; smallest-eigval eigenvector."""
    m = lset_map(alm, lset, lmax); w = m * m
    M = np.einsum('p,ip,jp->ij', w, VEC, VEC)
    ev, evec = np.linalg.eigh(M)
    n = evec[:, 0]
    return n / np.linalg.norm(n)

def ang(a, b):
    return float(np.degrees(np.arccos(np.clip(abs(np.dot(a, b)), 0, 1))))


# ------------------------------------------------------------------ five-fold ledger statistic
def fivefold_H(mp, nhat, nb=16):
    """Convention-free fractional m=5 azimuthal power about nhat, ring-averaged in polar angle."""
    nhat = np.asarray(nhat, float); nhat = nhat / np.linalg.norm(nhat)
    a = np.array([0., 0., 1.]) if abs(nhat[2]) < 0.9 else np.array([1., 0., 0.])
    e1 = a - (a @ nhat) * nhat; e1 /= np.linalg.norm(e1); e2 = np.cross(nhat, e1)
    pz = VEC.T @ nhat
    psi = np.arctan2(VEC.T @ e2, VEC.T @ e1)
    beta = np.arccos(np.clip(pz, -1, 1))
    edges = np.linspace(0, np.pi, nb + 1); H = 0.0; ct = 0
    for i in range(nb):
        sel = (beta >= edges[i]) & (beta < edges[i + 1])
        if sel.sum() < 20:
            continue
        f = mp[sel]
        c5 = np.mean(f * np.exp(-1j * 5 * psi[sel]))
        p0 = np.mean(f * f) + 1e-18
        H += abs(c5) ** 2 / p0; ct += 1
    return H / max(ct, 1)


# ------------------------------------------------------------------ parity statistic
def parity_g_from_cl(clv, lmax):
    """g = P+/P- = even-l power / odd-l power in D_l, 2<=l<=lmax.  g<1 => odd-parity preference."""
    ll = np.arange(2, lmax + 1)
    D = ll * (ll + 1) * clv[2:lmax + 1] / (2 * np.pi)
    ev = D[(ll % 2 == 0)].sum(); od = D[(ll % 2 == 1)].sum()
    return ev / od

def parity_g_from_alm(alm, lmax_stat, lmax=LMAX):
    cl = hp.alm2cl(alm, lmax=lmax)
    return parity_g_from_cl(cl, lmax_stat)


# ================================================================== THE PIPELINE
def gauss_beam(lmax, fwhm_deg):
    return hp.gauss_beam(np.radians(fwhm_deg), lmax=lmax)

def make_mask(nside=NSIDE, gal_cut_deg=20.0, apod_deg=5.0, seed_holes=0, rng=None):
    """Synthetic analysis mask: symmetric Galactic-plane cut, cosine-apodized.
       Replace with a real Planck common mask load for a data run."""
    theta, phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    b = 90.0 - np.degrees(theta)                       # 'latitude'
    mask = np.ones(hp.nside2npix(nside))
    inside = np.abs(b) < gal_cut_deg
    mask[inside] = 0.0
    if apod_deg > 0:                                   # cosine taper across apod band
        band = (np.abs(b) >= gal_cut_deg) & (np.abs(b) < gal_cut_deg + apod_deg)
        x = (np.abs(b[band]) - gal_cut_deg) / apod_deg
        mask[band] = 0.5 * (1 - np.cos(np.pi * x))
    return mask

def _binary(mask):
    return (mask > 0.5).astype(float)

def inpaint(mp, mask, lmax=LMAX, niter=40):
    """Iterative harmonic ('cooling') inpainting. Applied IDENTICALLY to null and data."""
    m = _binary(mask)
    filled = mp * m
    for _ in range(niter):
        alm = hp.map2alm(filled, lmax=lmax, iter=0)
        rec = hp.alm2map(alm, NSIDE)
        filled = m * mp + (1 - m) * rec              # keep data outside mask, reconstruction inside
    return filled

def observe(alm_true, fwhm_deg=5.0, mask=None, noise_uK=1.0, inpaint_iter=40,
            rng=None, lmax=LMAX):
    """
    Carry a TRUE alm through the full chain and return (processed_map, recovered_alm).
    Steps: beam -> map -> homogeneous white noise -> apply mask -> inpaint -> recover alm.
    The SAME function is used for null sims and for data, guaranteeing an identical pipeline.
    """
    if rng is None:
        rng = np.random.default_rng()
    bl = gauss_beam(lmax, fwhm_deg)
    alm_b = hp.almxfl(alm_true.copy(), bl)
    mp = hp.alm2map(alm_b, NSIDE)
    if noise_uK and noise_uK > 0:
        mp = mp + rng.normal(0.0, noise_uK, size=mp.size)
    if mask is None:
        recov = hp.map2alm(mp, lmax=lmax, iter=1)
        return mp, recov
    mp_ip = inpaint(mp, mask, lmax=lmax, niter=inpaint_iter)
    recov = hp.map2alm(mp_ip, lmax=lmax, iter=1)
    return mp_ip, recov

def draw_true_alm(cl=None, lmax=LMAX, rng=None):
    if cl is None:
        cl = get_cl(lmax)
    # healpy synalm uses its own RNG; for reproducibility across a run we seed via numpy state
    if rng is not None:
        st = np.random.get_state()
        np.random.seed(int(rng.integers(0, 2**31 - 1)))
        alm = hp.synalm(cl, lmax=lmax)
        np.random.set_state(st)
        return alm
    return hp.synalm(cl, lmax=lmax)


# ------------------------------------------------------------------ p-value helpers
def one_sided_p(null_vals, obs, tail='low'):
    null_vals = np.asarray(null_vals, float)
    if tail == 'low':
        return float((null_vals <= obs).mean())
    return float((null_vals >= obs).mean())

def zscore(null_vals, obs):
    null_vals = np.asarray(null_vals, float)
    return float((obs - null_vals.mean()) / (null_vals.std() + 1e-30))


if __name__ == "__main__":
    cl = get_cl()
    print("eim_pipeline self-check")
    print("  Cl[2,3,5,10] =", [round(float(cl[i]), 2) for i in (2, 3, 5, 10)])
    rng = np.random.default_rng(0)
    mask = make_mask()
    fsky = _binary(mask).mean()
    print("  mask f_sky = %.3f (gal cut 20 deg)" % fsky)
    a = draw_true_alm(cl, rng=rng)
    mp, rec = observe(a, mask=mask, rng=rng)
    print("  observe() ok: map std=%.2f  recovered alm size=%d" % (mp.std(), rec.size))
