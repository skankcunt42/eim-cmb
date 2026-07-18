"""
parity_test.py  --  g_L = P+/P- parity ratio, HARDENED.

WHAT CHANGED vs the original:
  (1) MATCHED NULL. Original null was a bare per-l chi^2 on full-sky theory C_l. Replaced by the
      shared pipeline (beam + mask + homogeneous noise + inpainting): each null sim is carried
      through eim_pipeline.observe() and g_L is read from its recovered pseudo-C_l. Mask coupling,
      which biases even/odd power, is now inside the null.
  (2) LOOK-ELSEWHERE MADE EXPLICIT AND CORRECTED. The predeclaration forbids scanning L to minimise
      p. This script (a) reports the LOCAL p at the single PREDECLARED L0, and (b) separately
      quantifies the look-elsewhere inflation by calibrating the distribution of min-over-L local p
      on the null, giving a GLOBAL (scan-corrected) p. The scan number is shown only to expose the
      inflation, never as the headline.

DECISION NUMBER = local p at the predeclared L0. Global p is the honest cost if one had scanned.
"""
import warnings; warnings.simplefilter('ignore')
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import healpy as hp
import eim_pipeline as E
from eim_pipeline import get_cl, draw_true_alm, observe, make_mask, parity_g_from_cl

rng = np.random.default_rng(2026)
LMAX_PIPE = 48                      # analysis band for parity
L0 = 20                            # THE predeclared decision L (fixed before data; do not move)
LSCAN = np.arange(6, 41, 2)        # scan range, used ONLY to measure look-elsewhere inflation
FWHM, NOISE = 5.0, 1.0
mask = make_mask(gal_cut_deg=20.0, apod_deg=5.0)
NMC = 800
cl = get_cl(LMAX_PIPE)

def g_at(clv, L): return parity_g_from_cl(clv, L)

# ---------------- MATCHED NULL: pipeline -> recovered pseudo-C_l -> g over all L ----------------
gsim = np.zeros((NMC, len(LSCAN)))
for i in range(NMC):
    a = draw_true_alm(cl, lmax=LMAX_PIPE, rng=rng)
    _, rec = observe(a, fwhm_deg=FWHM, mask=mask, noise_uK=NOISE, rng=rng, lmax=LMAX_PIPE)
    clr = hp.alm2cl(rec, lmax=LMAX_PIPE)
    for j, L in enumerate(LSCAN):
        gsim[i, j] = g_at(clr, L)
j0 = int(np.where(LSCAN == L0)[0][0])
mean, std = gsim.mean(0), gsim.std(0)
g_theory = parity_g_from_cl(cl, L0)
print("REFERENCE FLAW FIXED: the LCDM-shape expectation of this g estimator is %.3f, NOT 1.0" % g_theory)
print("   (even-l sum includes the large l=2 Sachs-Wolfe term; g=1 is the wrong null center).")
print("   Any 'odd-parity deficit' must be measured from the matched null below, not from g=1.\n")
print("MATCHED NULL g_L (through beam+mask+noise+inpaint), NMC=%d" % NMC)
print("   at predeclared L0=%d:  mean=%.3f  std=%.3f  (full-sky chi^2 null would give a different, tighter std)" %
      (L0, mean[j0], std[j0]))

# local p at each L: one-sided (odd-parity excess => low g)
def local_p_row(grow):
    # grow: shape (len(LSCAN),) observed g at each L; compare to null columns
    return np.array([(gsim[:, j] <= grow[j]).mean() for j in range(len(LSCAN))])

# ---------------- look-elsewhere calibration: distribution of min-over-L local p under the null ----
# For each null sim, compute its local p at each L against the OTHER sims, take the min => scan stat.
# Approx (leave-one-out is expensive): rank-based local p using the full null ensemble.
ranks = np.argsort(np.argsort(gsim, axis=0), axis=0)         # 0..NMC-1 per column
local_p_null = (ranks + 1) / (NMC + 1)                        # one-sided low-tail p per sim per L
minp_null = local_p_null.min(axis=1)                          # scan statistic under the null
def global_p(min_local_p): return float((minp_null <= min_local_p).mean())

print("   look-elsewhere: scanning L in [%d..%d] (%d points)." % (LSCAN[0], LSCAN[-1], len(LSCAN)))
print("   a LOCAL p=0.05 found by scanning corresponds to GLOBAL p=%.3f (inflation x%.1f)" %
      (global_p(0.05), global_p(0.05) / 0.05))

# ---------------- POWER: inject even-l suppression, report LOCAL@L0 and GLOBAL(scan) ----------------
def observed_g_row(A):
    """inject even-l suppression by A, push through pipeline, return g at every L (median of reps)."""
    reps = []
    for _ in range(12):
        ci = cl.copy(); ll = np.arange(len(ci));
        a = draw_true_alm(ci, lmax=LMAX_PIPE, rng=rng)
        # suppress even-l a_lm amplitude by sqrt(1-A) so C_l(even) *= (1-A)
        for L in range(2, LMAX_PIPE + 1):
            if L % 2 == 0:
                for m in range(L + 1):
                    a[hp.Alm.getidx(LMAX_PIPE, L, m)] *= np.sqrt(max(1 - A, 0))
        _, rec = observe(a, fwhm_deg=FWHM, mask=mask, noise_uK=NOISE, rng=rng, lmax=LMAX_PIPE)
        clr = hp.alm2cl(rec, lmax=LMAX_PIPE)
        reps.append([g_at(clr, L) for L in LSCAN])
    return np.median(np.array(reps), axis=0)

print("\nPOWER (inject even-l suppression A; report LOCAL p at L0=%d and GLOBAL scan-corrected p):" % L0)
for A in [0.0, 0.10, 0.20, 0.30]:
    grow = observed_g_row(A)
    lp = local_p_row(grow)
    p_local0 = lp[j0]
    p_glob = global_p(lp.min())
    z0 = (grow[j0] - mean[j0]) / std[j0]
    print("   A=%2.0f%%  g(L0)=%.3f  z(L0)=%+.2f  p_local@L0=%.4f   [scan min p_local=%.4f -> p_global=%.4f]" %
          (A * 100, grow[j0], z0, p_local0, lp.min(), p_glob))

# demo 'observed-like' odd preference value at L0 (real Planck pseudo-C_l plugs in here)
G_OBS_DEMO = 0.80
p_demo_local = float((gsim[:, j0] <= G_OBS_DEMO).mean())
z_demo = (G_OBS_DEMO - mean[j0]) / std[j0]
print("\nDEMO observed-like g(L0=%d)=%.2f -> z=%+.2f  p_local=%.4f (matched null; real Planck plugs in here)" %
      (L0, G_OBS_DEMO, z_demo, p_demo_local))

# ---------------- FIGURE ----------------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
ax[0].hist(gsim[:, j0], bins=50, color='#4C78A8', alpha=.85, density=True, label='matched ΛCDM null g(L0)')
ax[0].axvline(1.0, color='k', ls='--', lw=1, label='g=1')
ax[0].axvline(G_OBS_DEMO, color='#E45756', lw=2, label='demo observed-like g=0.80')
ax[0].set_xlabel('parity ratio g = P⁺/P⁻ at predeclared L0=%d' % L0); ax[0].set_ylabel('density')
ax[0].set_title('Parity null — beam+mask+noise+inpaint'); ax[0].legend(fontsize=8)
# look-elsewhere: local vs global p across the scan for the demo value applied at every L
lp_demo = np.array([(gsim[:, j] <= G_OBS_DEMO).mean() for j in range(len(LSCAN))])
ax[1].plot(LSCAN, lp_demo, '-o', color='#E45756', ms=4, label='local p (per L) for g=0.80')
ax[1].axhline(0.05, color='grey', ls=':', label='p=0.05')
ax[1].axvline(L0, color='k', ls='--', lw=1, label='predeclared L0')
ax[1].set_xlabel('ℓ_max (scan)'); ax[1].set_ylabel('local one-sided p'); ax[1].set_yscale('log')
ax[1].set_title('Look-elsewhere: p slides with L (why L0 is frozen)'); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('parity_result.png', dpi=95)
print("saved parity_result.png")
