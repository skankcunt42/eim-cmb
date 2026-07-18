"""
reproduce_paper_numbers.py -- referee-grade regeneration of the section-3 / Table-1 numbers of the
closure paper directly from the engine, with Monte-Carlo error bars and pass/flag tolerances.

Each row recomputes a paper number from eim_pipeline primitives through the matched pipeline
(beam + apodized symmetric Galactic mask + homogeneous noise + iterative-harmonic inpaint) so the
comparison is like-for-like with the manuscript. Numbers are stochastic; we print the MC mean +-
s.e. and mark PASS when the paper value lies within ~2 s.e. or a stated absolute tolerance.

Run:  python reproduce_paper_numbers.py [N]     (default N=400 sims; use 800+ to tighten)
"""
import warnings; warnings.simplefilter("ignore")
import sys, numpy as np, healpy as hp
import eim_pipeline as E
from eim_pipeline import (NSIDE, LMAX, get_cl, draw_true_alm, make_mask, observe, axis_of, ang,
                          parity_g_from_alm, parity_g_from_cl, fivefold_H, band_map, _binary)

def _report(name, val, paper, tol, se=None, extra=""):
    hit = abs(val - paper) <= tol
    se_s = (" +- %.3f" % se) if se is not None else ""
    print("   %-46s = %7.3f%s   paper %-7s  [%s]%s"
          % (name, val, se_s, ("%.3f" % paper), "PASS" if hit else "FLAG", extra))
    return hit

def main(N=400):
    rng = np.random.default_rng(20260718)
    cl = get_cl(LMAX)
    mask = make_mask()                              # apodized symmetric Galactic cut
    zpole = np.array([0, 0, 1.0])                   # mask symmetry axis ('Galactic pole' here)
    print("="*80)
    print("SECTION-3 NUMBER REPRODUCTION  (N=%d sims, matched pipeline, f_sky=%.3f)"
          % (N, _binary(mask).mean()))
    print("-"*80)

    # ---- collect per-sim statistics through the pipeline ----
    g_full, g_mask = [], []
    alpha_mask, alpha_full = [], []
    even_axis_to_pole = []
    H_null = []
    for _ in range(N):
        a = draw_true_alm(cl, rng=rng)
        # full-sky (no mask) recovery
        _, af = observe(a, mask=None, rng=rng)
        g_full.append(parity_g_from_alm(af, 20))
        ne_f, no_f = axis_of(af, [2]), axis_of(af, [3])
        alpha_full.append(ang(ne_f, no_f))
        # matched masked recovery
        _, am = observe(a, mask=mask, rng=rng)
        g_mask.append(parity_g_from_alm(am, 20))
        ne, no = axis_of(am, [2]), axis_of(am, [3])
        nJ = axis_of(am, [2, 3])
        alpha_mask.append(ang(ne, no))
        even_axis_to_pole.append(ang(ne, zpole))
        H_null.append(fivefold_H(band_map(am, 5, 12), nJ))

    g_full = np.array(g_full); g_mask = np.array(g_mask)
    alpha_mask = np.array(alpha_mask); alpha_full = np.array(alpha_full)
    even_axis_to_pole = np.array(even_axis_to_pole); H_null = np.array(H_null)
    def se_of_std(x): return x.std() / np.sqrt(2*len(x))    # approx s.e. of a std
    def se_of_med(x): return 1.253*x.std()/np.sqrt(len(x))  # approx s.e. of a median

    ok = []
    print("[3.1] Parity g_L = P+/P-")
    ok.append(_report("LCDM-shape reference g(L=20)", parity_g_from_cl(cl, 20), 1.126, 0.01))
    ok.append(_report("full-sky std of g(20)", g_full.std(), 0.198, 0.03, se_of_std(g_full)))
    ok.append(_report("masked std of g(20)", g_mask.std(), 0.27, 0.04, se_of_std(g_mask),
                      extra=" infl=%.0f%%" % (100*(g_mask.std()/g_full.std()-1))))
    ok.append(_report("matched-null mean g(20) (paper 1.17-1.20)", g_mask.mean(), 1.185, 0.04,
                      g_mask.std()/np.sqrt(N)))

    print("[3.2] Co-axiality alpha = arccos|n2 . n3|")
    ok.append(_report("matched-null median alpha (deg)", np.median(alpha_mask), 54.1, 3.0,
                      se_of_med(alpha_mask)))
    ok.append(_report("full-sky control median alpha (deg)", np.median(alpha_full), 59.93, 3.0,
                      se_of_med(alpha_full)))
    ok.append(_report("fraction alpha < 10 deg (matched)", (alpha_mask < 10).mean(), 0.020, 0.015))
    # frame audit: recovered even-axis angle to mask pole
    ok.append(_report("frame audit: median even-axis-to-pole (deg)", np.median(even_axis_to_pole),
                      73.6, 5.0, se_of_med(even_axis_to_pole)))
    # KS of cos(angle) against uniform (isotropy in cos)
    c = np.abs(np.cos(np.radians(even_axis_to_pole)))
    cs = np.sort(c); Femp = np.arange(1, len(cs)+1)/len(cs); Fan = cs   # uniform on [0,1] in |cos|
    ks_frame = float(np.max(np.abs(Femp - Fan)))
    ok.append(_report("frame audit: KS(|cos| vs uniform)", ks_frame, 0.26, 0.10))

    print("[3.3] Five-fold ledger H5 (matched null, about frozen axis)")
    ok.append(_report("H5 null median", np.median(H_null), 0.047, 0.015))
    ok.append(_report("H5 null 95th percentile", np.percentile(H_null, 95), 0.078, 0.025))

    # ---- H5 null-vs-null calibration: split the null, KS + false-positive rate ----
    half = N // 2
    ref, test = H_null[:half], H_null[half:]
    # one-sided matched-null p for each 'test' draw against 'ref'
    ps = np.array([(ref >= h).mean() for h in test])
    cs2 = np.sort(ps); Femp2 = np.arange(1, len(cs2)+1)/len(cs2)
    ks_p = float(np.max(np.abs(Femp2 - cs2)))
    fpr = float((ps <= 0.05).mean())
    print("[calibration] H5 null-vs-null p-values")
    print("   KS(p vs uniform) D = %.3f   (paper D=0.07)   |   FPR@0.05 = %.3f  (paper 0.04)"
          % (ks_p, fpr))

    print("-"*80)
    npass = sum(ok)
    print("SUMMARY: %d/%d section-3 numbers reproduced within tolerance (N=%d; increase N to tighten)."
          % (npass, len(ok), N))
    print("="*80)
    return npass, len(ok)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    main(N)
