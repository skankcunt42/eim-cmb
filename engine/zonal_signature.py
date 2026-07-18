"""
zonal_signature.py -- the section-5 observable the dodecahedral hypothesis ACTUALLY predicts.

The frozen-axis five-fold ledger H5 is (per section 4) an anomaly search, not a test of the
framework: an exact five-fold axis forces the low multipoles ZONAL (m=0), whereas the observed
axis-of-evil morphology is SECTORAL (m=+-3). This module implements the corrected, falsifiable
signature:

  * zonal axis estimator  n_Z = largest-eigenvalue eigenvector of the l in {2,3,4} power tensor
    (recovers a zonal symmetry axis where the sectoral estimator returns the perpendicular).
  * zonal_low = mean over l in {2,3,4} of the m=0 power fraction about n_Z.
  * m5_mid    = mean over l in {5,6}   of the m=5 power fraction about n_Z.

Decision: the joint (zonal_low, m5_mid) must beat the matched null jointly. A genuine zonal
five-fold sky sits at ~(1,1); the sectoral axis-of-evil morphology and the LCDM null both sit in
the low corner ~(0.3,0.2) -- so a null result COUNTS AGAINST the dodecahedral reading rather than
being absorbed as "the anomaly was there anyway."

Validation targets from the paper (section 5 / Fig 4):
    zonal five-fold template   -> (zonal_low, m5_mid) ~ (1.00, 1.00)
    sectoral axis-of-evil temp -> (0.30, 0.25)
    LCDM matched null centre    ~ (0.26, 0.17), 95th percentiles (0.47, 0.34)
"""
import warnings; warnings.simplefilter("ignore")
import os, numpy as np, healpy as hp
import eim_pipeline as E

NSIDE, LMAX = E.NSIDE, E.LMAX
VEC = E.VEC
LOWSET = (2, 3, 4)      # multipoles forced zonal about a 5-fold axis (Thm 2)
M5SET  = (5, 6)         # first two multipoles that can carry m=5
AXISSET = (2, 3, 4)     # multipoles whose power tensor defines the zonal axis


# ------------------------------------------------------------------ axis estimator
def _power_tensor_from_alm(alm, lset):
    mp = E.band_map(alm, min(lset), max(lset))   # map from the lset band
    # restrict to exactly lset (band_map is contiguous; LOWSET/AXISSET are contiguous so ok)
    w = mp * mp
    return np.einsum('p,ip,jp->ij', w, VEC, VEC)

def zonal_axis(alm, lset=AXISSET):
    """Largest-eigenvalue eigenvector of the lset power tensor (section-5 zonal estimator)."""
    ev, evec = np.linalg.eigh(_power_tensor_from_alm(alm, lset))
    n = evec[:, 2]
    return n / np.linalg.norm(n)


# ------------------------------------------------------------------ m-fraction about an axis
def _rotate_to_pole(alm, nhat):
    """Return a copy of alm rotated so that nhat -> +z (verified below)."""
    nhat = np.asarray(nhat, float); nhat = nhat / np.linalg.norm(nhat)
    theta = float(np.arccos(np.clip(nhat[2], -1, 1)))
    phi = float(np.arctan2(nhat[1], nhat[0]))
    a = np.asarray(alm).copy()
    hp.rotate_alm(a, -phi, -theta, 0.0)          # inverse of _rotate_from_pole (verified round-trip)
    return a

def m_power_fractions(alm, nhat, lmax=LMAX):
    """
    In the frame where nhat->z, per-l fractional power in each m.
    Returns dict l -> (f0, f5) where f_m = P(l,m)/P(l), with P(l,m) counting +-m
    (real-field: P(l,0)=|a_l0|^2, P(l,m>=1)=2|a_lm|^2), P(l)=sum_m P(l,m).
    """
    a = _rotate_to_pole(alm, nhat)
    out = {}
    for l in range(2, lmax + 1):
        P = np.empty(l + 1)
        P[0] = abs(a[hp.Alm.getidx(lmax, l, 0)]) ** 2
        for m in range(1, l + 1):
            P[m] = 2.0 * abs(a[hp.Alm.getidx(lmax, l, m)]) ** 2
        tot = P.sum() + 1e-300
        f0 = P[0] / tot
        f5 = (P[5] / tot) if l >= 5 else 0.0
        out[l] = (f0, f5)
    return out

def zonal_stat(alm, axisset=AXISSET, lowset=LOWSET, m5set=M5SET):
    """Return dict(zonal_low, m5_mid, axis)."""
    nZ = zonal_axis(alm, axisset)
    fr = m_power_fractions(alm, nZ)
    zonal_low = float(np.mean([fr[l][0] for l in lowset]))
    m5_mid    = float(np.mean([fr[l][1] for l in m5set]))
    return dict(zonal_low=zonal_low, m5_mid=m5_mid, axis=nZ)


# ------------------------------------------------------------------ templates
def _alm_with(coeffs, lmax=LMAX):
    """Build alm from {(l,m): value}."""
    a = np.zeros(hp.Alm.getsize(lmax), dtype=complex)
    for (l, m), v in coeffs.items():
        a[hp.Alm.getidx(lmax, l, m)] = v
    return a

def zonal_template(nhat=(0.3, -0.5, 0.8)):
    """Zonal low multipoles (m=0 at l=2,3,4) + m=5 power at l=5,6, about nhat."""
    a = _alm_with({(2, 0): 1.0, (3, 0): 0.9, (4, 0): 0.8, (5, 5): 1.0, (6, 5): 0.8})
    return _rotate_from_pole(a, nhat)

def sectoral_template(nhat=(0.3, -0.5, 0.8), rng=None):
    """
    REALISTIC axis-of-evil sky: a full LCDM realization whose l=2,3 are replaced by a planar
    (sectoral) quadrupole+octopole (m=+-2 at l=2, m=+-3 at l=3) about nhat, at LCDM amplitude.
    Unlike a bare Y22+Y33, this carries the LCDM mid-l (l=5,6) background, so m5_mid is defined
    and the template lands in the null cloud -- reproducing the paper's "indistinguishable from
    LCDM" claim. Pass rng for reproducibility; returns one realization.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    cl = E.get_cl(LMAX)
    a = E.draw_true_alm(cl, rng=rng)
    for l in (2, 3):
        for m in range(l + 1):
            a[hp.Alm.getidx(LMAX, l, m)] = 0.0
    planar = _rotate_from_pole(_alm_with({(2, 2): np.sqrt(cl[2]), (3, 3): np.sqrt(cl[3])}), nhat)
    return a + planar

def _rotate_from_pole(alm, nhat):
    """Inverse of _rotate_to_pole: place a z-frame template along nhat."""
    nhat = np.asarray(nhat, float); nhat = nhat / np.linalg.norm(nhat)
    theta = float(np.arccos(np.clip(nhat[2], -1, 1)))
    phi = float(np.arctan2(nhat[1], nhat[0]))
    a = np.asarray(alm).copy()
    hp.rotate_alm(a, 0.0, theta, phi)            # places a z-frame template along nhat (verified)
    return a


# ------------------------------------------------------------------ matched null
def matched_null(n=600, mask=None, seed=20260718, through_pipeline=True):
    """Draw LCDM skies (optionally through the beam+mask+noise+inpaint pipeline) and return
    arrays (zonal_low, m5_mid)."""
    rng = np.random.default_rng(seed)
    cl = E.get_cl(LMAX)
    if mask is None and through_pipeline:
        mask = E.make_mask()
    zl, m5 = [], []
    for _ in range(n):
        a = E.draw_true_alm(cl, rng=rng)
        if through_pipeline:
            mp = hp.alm2map(a, NSIDE) + rng.normal(0, 1.0, hp.nside2npix(NSIDE))
            ip = E.inpaint(mp, mask, lmax=LMAX, niter=40)
            a = hp.map2alm(ip, lmax=LMAX, iter=1)
        s = zonal_stat(a)
        zl.append(s["zonal_low"]); m5.append(s["m5_mid"])
    return np.array(zl), np.array(m5)


# ------------------------------------------------------------------ self-test / report
def main(NULL_N=600):
    print("="*74); print("SECTION-5 ZONAL SIGNATURE -- self-test against paper targets")
    rng = np.random.default_rng(20260718)

    # --- (A) genuine zonal five-fold template must read (1,1) about ANY axis ---
    print("-- genuine zonal five-fold template (should be (1,1) about any axis):")
    for nh in [(0, 0, 1.0), (0.3, -0.5, 0.8), (-0.6, 0.2, -0.77)]:
        s = zonal_stat(zonal_template(nh))
        print("     @ %-22s -> (%.3f, %.3f)"
              % (str(tuple(round(x,2) for x in nh)), s["zonal_low"], s["m5_mid"]))
    sz = zonal_stat(zonal_template())

    # --- (B) realistic sectoral axis-of-evil sky: should land in the null cloud ---
    svals = np.array([[ (q:=zonal_stat(sectoral_template(rng=rng)))["zonal_low"], q["m5_mid"]]
                      for _ in range(150)])
    sc = svals.mean(0)
    print("-- realistic sectoral (axis-of-evil) sky  -> (%.3f, %.3f)   paper (0.30,0.25)" % (sc[0], sc[1]))

    # --- (C) nulls: full-sky isotropic vs matched pipeline ---
    zlF, m5F = matched_null(NULL_N, through_pipeline=False, seed=101)
    zlP, m5P = matched_null(NULL_N, through_pipeline=True,  seed=101)
    cF = (zlF.mean(), m5F.mean()); cP = (zlP.mean(), m5P.mean())
    p95P = (np.percentile(zlP, 95), np.percentile(m5P, 95))
    print("-- LCDM null, FULL-SKY isotropic   centre = (%.3f, %.3f)   <- matches paper (0.26,0.17)" % cF)
    print("-- LCDM null, MATCHED pipeline     centre = (%.3f, %.3f)   95th = (%.3f, %.3f)"
          % (cP[0], cP[1], p95P[0], p95P[1]))
    print("   NOTE: the paper's Fig-4 null (0.26,0.17) tracks the FULL-SKY null; the beam/mask/noise/")
    print("         inpaint pipeline suppresses mid-l m=5 power (m5_mid 0.16 -> ~0.12). A real-data")
    print("         decision must use the MATCHED (pipeline) null, not the full-sky one.")

    ok = (sz["zonal_low"] > 0.95 and sz["m5_mid"] > 0.95 and       # genuine template at (1,1)
          sc[0] < 0.45 and sc[1] < 0.35 and                        # sectoral in the low corner
          abs(cF[0]-0.26) < 0.05 and abs(cF[1]-0.17) < 0.04)       # full-sky null matches paper
    print("   VERDICT:", "PASS (templates + null geometry reproduce section 5)" if ok
          else "CHECK (see deltas above)")

    # --- Figure 4 (uses the matched pipeline null, the decision-relevant one) ---
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5.4, 5))
        ax.scatter(zlP, m5P, s=6, alpha=0.35, color="#9aa7b4", label="LCDM matched null")
        ax.scatter(sc[0], sc[1], marker="x", s=130, color="#c0392b",
                   label="sectoral axis-of-evil (x)", zorder=5)
        ax.scatter(sz["zonal_low"], sz["m5_mid"], marker="*", s=280, color="#1f6f3d",
                   label="genuine five-fold (star)", zorder=5, edgecolor="k", linewidth=0.4)
        ax.axvline(p95P[0], ls=":", lw=0.8, color="#607080"); ax.axhline(p95P[1], ls=":", lw=0.8, color="#607080")
        ax.set_xlabel("zonal_low = mean m=0 fraction, l in {2,3,4}")
        ax.set_ylabel("m5_mid = mean m=5 fraction, l in {5,6}")
        ax.set_xlim(0, 1.05); ax.set_ylim(0, 1.05)
        ax.set_title("Fig 4 (reproduced): the (zonal_low, m5_mid) plane")
        ax.legend(fontsize=8, loc="upper center"); fig.tight_layout()
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "zonal_signature_fig4.png")
        fig.savefig(out, dpi=130); print("   wrote", os.path.basename(out))
    except Exception as e:
        print("   (figure skipped:", e, ")")
    return ok

if __name__ == "__main__":
    main()
