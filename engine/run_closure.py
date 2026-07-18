"""
run_closure.py -- full-closure driver for the frozen-axis protocol across >=2 component-separated
maps AND >=2 masks, including BOTH the (reclassified) five-fold ledger H5 and the section-5 zonal
signature. Extends run_real.py from a >=2-map first look to the decision rule of predeclaration
section 2.

DATA MODE (default): reads {smica,nilc,sevem,commander}_nside64_uK.fits + a real mask produced by
    fetch_planck_streamlined.ipynb / fetch_planck.py, if present in the working dir.
DRY-RUN MODE (auto): if no maps are found, synthesizes LCDM 'pseudo-data' so the driver, the
    per-mask matched nulls, and the decision logic are all exercised end-to-end without archive
    access. Clearly labelled as NOT a science result.

Masks: the pipeline needs >=2. We use (A) the real common mask if present else a 20-deg apodized
cut, and (B) a more aggressive 30-deg apodized cut. For a submission, replace mask B with a SECOND
REAL Planck mask (e.g. a confidence/UTA76 variant) -- add it to the fetch and drop it in here.

Decision rules:
  * legacy anomaly set {g_L, alpha, H5}: joint beat of the matched null across >=2 maps AND
    >=2 masks, look-elsewhere corrected, with H5 locked to a common frozen axis. Per the paper H5
    is an ANOMALY search, not a framework test.
  * section-5 framework test {zonal_low, m5_mid}: both statistics jointly beat the matched null
    (>= its 95th percentile) across >=2 maps AND >=2 masks, about a common ZONAL axis.
"""
import warnings; warnings.simplefilter("ignore")
import os, numpy as np, healpy as hp
import eim_pipeline as E
from eim_pipeline import (LMAX, NSIDE, get_cl, draw_true_alm, make_mask, axis_of, ang,
                          band_map, fivefold_H, parity_g_from_alm, parity_g_from_cl, inpaint,
                          observe, _binary)
import zonal_signature as Z

L0, BAND = 20, (5, 12)
MAP_FILES = [("SMICA", "smica_nside64_uK.fits"), ("NILC", "nilc_nside64_uK.fits"),
             ("SEVEM", "sevem_nside64_uK.fits"), ("COMMANDER", "commander_nside64_uK.fits")]


def _load_binary(path):
    m = hp.read_map(path, dtype=np.float64)
    if hp.get_nside(m) != NSIDE:
        m = hp.ud_grade(m, NSIDE)
    return (m > 0.5).astype(float)

def load_masks():
    """>=2 masks. Prefer TWO real Planck masks (common + GAL070) from the streamlined fetch;
    fall back to synthetic cuts (with a warning) for any that are missing."""
    masks = []
    # mask 1: real common mask, else synthetic 20deg
    if os.path.exists("commonmask_nside64.fits"):
        masks.append(("common", _load_binary("commonmask_nside64.fits")))
    else:
        print("  [masks] commonmask_nside64.fits not found -> synthetic 20deg cut (NOT a real mask)")
        masks.append(("synthGC20", make_mask(gal_cut_deg=20.0)))
    # mask 2: real GAL070 Galactic-plane mask, else synthetic 30deg
    if os.path.exists("gal070mask_nside64.fits"):
        masks.append(("gal070", _load_binary("gal070mask_nside64.fits")))
    else:
        print("  [masks] gal070mask_nside64.fits not found -> synthetic 30deg cut (NOT a real 2nd")
        print("          mask). Re-run fetch_planck_streamlined.ipynb to get the real GAL070 mask.")
        masks.append(("synthGC30", make_mask(gal_cut_deg=30.0)))
    return masks


def load_data(rng):
    """Return list of (name, alm_or_map) for available maps, or synthesize LCDM pseudo-data."""
    have = [(n, f) for n, f in MAP_FILES if os.path.exists(f)]
    if have:
        out = []
        for n, f in have:
            m = hp.read_map(f, dtype=np.float64)
            if hp.get_nside(m) != NSIDE: m = hp.ud_grade(m, NSIDE)
            out.append((n, m))
        return out, True
    # dry-run: synthesize 4 correlated LCDM skies (shared low-l + independent noise) as pseudo-data
    cl = get_cl(LMAX)
    base = draw_true_alm(cl, rng=rng)
    out = []
    for n, _ in MAP_FILES:
        a = base.copy() + 0.15 * draw_true_alm(cl, rng=rng)      # near-common sky, small differences
        out.append((n, hp.alm2map(a, NSIDE)))
    return out, False


def matched_null_all(mask, n, rng):
    """Matched null through `mask` for all statistics. Returns dict of arrays."""
    cl = get_cl(LMAX)
    g, al, H, zl, m5 = [], [], [], [], []
    for _ in range(n):
        a = draw_true_alm(cl, rng=rng)
        mp = hp.alm2map(a, NSIDE) + rng.normal(0, 1.0, hp.nside2npix(NSIDE))
        ip = inpaint(mp, mask, lmax=LMAX, niter=40)
        rec = hp.map2alm(ip, lmax=LMAX, iter=1)
        nJ = axis_of(rec, [2, 3])
        g.append(parity_g_from_alm(rec, L0))
        al.append(ang(axis_of(rec, [2]), axis_of(rec, [3])))
        H.append(fivefold_H(band_map(rec, *BAND), nJ))
        s = Z.zonal_stat(rec)
        zl.append(s["zonal_low"]); m5.append(s["m5_mid"])
    return {k: np.array(v) for k, v in dict(g=g, alpha=al, H=H, zl=zl, m5=m5).items()}


def data_stats(mp, mask):
    ip = inpaint(mp, mask, lmax=LMAX, niter=40)
    rec = hp.map2alm(ip, lmax=LMAX, iter=1)
    nJ = axis_of(rec, [2, 3])
    s = Z.zonal_stat(rec)
    return dict(g=parity_g_from_alm(rec, L0),
                alpha=ang(axis_of(rec, [2]), axis_of(rec, [3])),
                H=fivefold_H(band_map(rec, *BAND), nJ),
                nJ=nJ, zl=s["zonal_low"], m5=s["m5_mid"], nZ=s["axis"])


def gal_lb(n):
    th, ph = hp.vec2ang(np.asarray(n, float))
    th = float(np.atleast_1d(th)[0]); ph = float(np.atleast_1d(ph)[0])
    return float(np.degrees(ph)), float(90.0 - np.degrees(th))


def main(NULL_N=300):
    rng = np.random.default_rng(20260718)
    data, is_real = load_data(rng)
    masks = load_masks()
    print("="*82)
    print("EIM-CMB FULL-CLOSURE RUN  |  %s  |  maps=%d  masks=%d  null N=%d"
          % ("REAL DATA" if is_real else "DRY-RUN (LCDM pseudo-data, NOT a science result)",
             len(data), len(masks), NULL_N))
    if not is_real:
        print("  (no *_nside64_uK.fits found; run fetch_planck_streamlined.ipynb first for real maps)")
    print("="*82)

    # results[(map,mask)] = dict of stats+p-values ; axes for cross-map stability
    results, axes_J, axes_Z = {}, {}, {}
    for mname, mask in masks:
        null = matched_null_all(mask, NULL_N, rng)
        p95_zl, p95_m5 = np.percentile(null["zl"], 95), np.percentile(null["m5"], 95)
        print("\n### MASK = %s  (f_sky=%.3f) ###" % (mname, _binary(mask).mean()))
        print("    null: g mean=%.3f | alpha med=%.1f | H5 med=%.4f 95th=%.4f | "
              "zonal_low 95th=%.3f | m5_mid 95th=%.3f"
              % (null["g"].mean(), np.median(null["alpha"]), np.median(null["H"]),
                 np.percentile(null["H"], 95), p95_zl, p95_m5))
        for dname, mp in data:
            s = data_stats(mp, mask)
            p_g = float((null["g"] <= s["g"]).mean())        # odd-parity excess -> low g
            p_a = float((null["alpha"] <= s["alpha"]).mean())  # co-axial -> small alpha
            p_H = float((null["H"] >= s["H"]).mean())          # more five-fold -> high H
            p_zl = float((null["zl"] >= s["zl"]).mean())       # zonal -> high zonal_low
            p_m5 = float((null["m5"] >= s["m5"]).mean())       # five-fold -> high m5_mid
            results[(dname, mname)] = dict(**s, p_g=p_g, p_a=p_a, p_H=p_H, p_zl=p_zl, p_m5=p_m5)
            axes_J[(dname, mname)] = s["nJ"]; axes_Z[(dname, mname)] = s["nZ"]
            lJb, bJb = gal_lb(s["nJ"])
            print("  %-9s | g=%.3f p=%.3f | a=%4.1f p=%.3f | H5=%.4f p=%.3f || "
                  "zonal_low=%.3f p=%.3f  m5_mid=%.3f p=%.3f  nZ(l,b)=(%.0f,%+.0f)"
                  % (dname, s["g"], p_g, s["alpha"], p_a, s["H"], p_H,
                     s["zl"], p_zl, s["m5"], p_m5, *gal_lb(s["nZ"])))

    # ---------- decision aggregation ----------
    def joint_hits_legacy(alpha_p=0.05):
        # a (map,mask) 'hit' = all three of g,alpha,H beat null at alpha_p
        return {k for k, r in results.items()
                if r["p_g"] < alpha_p and r["p_a"] < alpha_p and r["p_H"] < alpha_p}

    def joint_hits_zonal(alpha_p=0.05):
        return {k for k, r in results.items() if r["p_zl"] < alpha_p and r["p_m5"] < alpha_p}

    def spans(hits):
        maps = {k[0] for k in hits}; msk = {k[1] for k in hits}
        return len(maps), len(msk)

    print("\n" + "="*82)
    hitsL = joint_hits_legacy(); mL, kL = spans(hitsL)
    hitsZ = joint_hits_zonal();  mZ, kZ = spans(hitsZ)
    # cross-map axis agreement (max pairwise separation of frozen axes among legacy hits)
    def max_axis_spread(hits, axmap):
        hs = list(hits)
        if len(hs) < 2: return None
        return max(ang(axmap[a], axmap[b]) for i, a in enumerate(hs) for b in hs[i+1:])
    sprJ = max_axis_spread(hitsL, axes_J); sprZ = max_axis_spread(hitsZ, axes_Z)

    print("LEGACY anomaly set {g, alpha, H5}  (H5 is an anomaly search, not a framework test):")
    print("   joint hits (p<0.05) span %d maps x %d masks%s"
          % (mL, kL, "" if sprJ is None else "  | frozen-axis spread <=%.1f deg" % sprJ))
    passL = (mL >= 2 and kL >= 2 and (sprJ is not None and sprJ < 10))
    print("   closure rule (>=2 maps AND >=2 masks, common axis):",
          "MET" if passL else "NOT met")

    print("SECTION-5 framework test {zonal_low, m5_mid} about the zonal axis:")
    print("   joint hits (p<0.05) span %d maps x %d masks%s"
          % (mZ, kZ, "" if sprZ is None else "  | zonal-axis spread <=%.1f deg" % sprZ))
    passZ = (mZ >= 2 and kZ >= 2 and (sprZ is not None and sprZ < 10))
    print("   framework rule (>=2 maps AND >=2 masks, common zonal axis):",
          "MET -> zonal five-fold signature detected" if passZ
          else "NOT met -> consistent with LCDM / no zonal five-fold signature")
    print("="*82)
    print("NOTE: apply look-elsewhere correction to the parity L-scan (see parity_test.py) and, for")
    print("      a submission, swap mask B for a second REAL Planck mask and the iterative inpainter")
    print("      for the survey map-maker's inpainter. H5 rows are reported as anomaly documentation.")
    return dict(passL=passL, passZ=passZ, is_real=is_real)

if __name__ == "__main__":
    import sys
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 300)
