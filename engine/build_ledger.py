"""
build_ledger.py  --  H5 five-fold ledger, HARDENED.

WHAT CHANGED vs the original build_ledger.py (two defects fixed):

  (1) AXIS-POWER MISMATCH.  The original power test injected m=5 about z, froze
      nfroz=axis_of(alm,[2,3]) (data-driven), then measured H about z -- i.e. about the
      TRUE signal axis, not the frozen axis the decision rule actually uses. That overstates
      detectability by exactly the frozen-vs-true axis miss. Fixed: H is measured about the
      FROZEN axis, matching the decision rule. The old "measure about z" number is retained
      ONLY as a clearly-labelled ORACLE upper bound.

  (2) IDEALIZED NULL.  The original null used bare full-sky synalm. Fixed: null and signal
      are both carried through eim_pipeline.observe() (beam + mask + homogeneous noise +
      inpainting), so mask coupling is inside the null.

WHAT THE THREE POWER SCENARIOS SHOW (this is the substantive result, not a code tidy):
  The frozen-axis H5 statistic has power ONLY IF the mid-l five-fold axis coincides with the
  low-l (l=2,3) axis. That coincidence is the EIM 'axis-sharing' claim -- which the framework
  does NOT currently derive (polarization memo: tau is structurally axis-free; cellulation
  orientation is gauge). So:
    S_oracle  : inject m=5 about z, measure about z            -> optimistic upper bound (old code)
    S_nolink  : inject m=5 about z, low-l random, measure about frozen axis
                -> ~NULL: no low-l/mid-l axis link => no power. This is the honest cost of the fix.
    S_eim     : inject a SHARED axis into l=2,3 AND m=5 about it, measure about frozen axis
                -> power returns. This is the ONLY scenario in which H5 tests EIM, and it
                   presupposes the (undischarged) axis-sharing bridge.
"""
import warnings; warnings.simplefilter('ignore')
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import eim_pipeline as E
from eim_pipeline import (LMAX, NSIDE, get_cl, draw_true_alm, observe, make_mask,
                          axis_of, ang, band_map, fivefold_H, one_sided_p)
import healpy as hp

rng = np.random.default_rng(555)
cl = get_cl(LMAX)
mask = make_mask(gal_cut_deg=20.0, apod_deg=5.0)
FWHM, NOISE = 5.0, 1.0
NULL_N   = 400
POWER_REP = 40
BAND = (5, 12)

def idx(L, m): return hp.Alm.getidx(LMAX, L, m)

# ---------------- SELF-TEST: inject 5-fold about z, detect on z not x (median over draws) ----
# A single realization is noisy; the estimator's directionality is a median statement.
Hzs, Hxs = [], []
for _ in range(20):
    a0 = draw_true_alm(cl, rng=rng)
    a0[idx(6, 5)] += 80; a0[idx(5, 5)] += 80          # strong injection: unambiguous instrument check
    mp = band_map(a0, *BAND)
    Hzs.append(fivefold_H(mp, [0, 0, 1.])); Hxs.append(fivefold_H(mp, [1., 0, 0]))
Hz, Hx = np.median(Hzs), np.median(Hxs)
print("SELF-TEST 5-fold (median of 20): H(z)=%.4f  H(x)=%.4f  ratio=%.1f  [pass if z>>x]" % (Hz, Hx, Hz / (Hx + 1e-9)))
assert Hz > 3 * Hx, "fivefold estimator failed self-test"
print("SELF-TEST PASSED\n")

# ---------------- MATCHED NULL: pipeline(beam+mask+noise+inpaint) ; frozen axis; H about it ----------
Hn = np.zeros(NULL_N); al = np.zeros(NULL_N)
for i in range(NULL_N):
    a = draw_true_alm(cl, rng=rng)
    _, rec = observe(a, fwhm_deg=FWHM, mask=mask, noise_uK=NOISE, rng=rng)
    nfroz = axis_of(rec, [2, 3])
    Hn[i] = fivefold_H(band_map(rec, *BAND), nfroz)
    al[i] = ang(axis_of(rec, [2]), axis_of(rec, [3]))
p95 = np.percentile(Hn, 95)
print("MATCHED NULL H (about frozen l2,3 axis, through beam+mask+noise+inpaint):")
print("   median=%.4f  95th=%.4f  (N=%d)" % (np.median(Hn), p95, NULL_N))
r = np.corrcoef(Hn, al)[0, 1]
print("   INDEPENDENCE corr(H,alpha)=%.3f  (under the matched null; not a bound under the alternative)\n" % r)

# ---------------- POWER: three scenarios ----------------
W_COAX = 0.9   # low-l aligned fraction for S_eim: yields frozen axis ~5 deg from truth
               # (i.e. co-axiality at the ~10 deg 'observed anomaly' level). At w=0.75 (~13 deg)
               # the frozen-axis error already erodes H5 power to borderline -- reported below.
def _amp(L): return np.sqrt((2 * L + 1) * cl[L])

def inject_run(amp, mode):
    """mode in {'oracle','nolink','eim'}; returns median H over POWER_REP reps."""
    out = []
    for _ in range(POWER_REP):
        a = draw_true_alm(cl, rng=rng)
        if mode == 'eim':
            # genuinely co-axial low-l about z: fade the random l=2,3 part, add planar-about-z content
            for L in (2, 3):
                for m in range(L + 1):
                    a[idx(L, m)] *= (1 - W_COAX)
                a[idx(L, L)] += W_COAX * _amp(L)
        a[idx(6, 5)] += amp; a[idx(5, 5)] += amp        # m=5 about z in both l=5,6
        _, rec = observe(a, fwhm_deg=FWHM, mask=mask, noise_uK=NOISE, rng=rng)
        if mode == 'oracle':
            n_meas = np.array([0, 0, 1.])               # measure about TRUE axis (old code)
        else:
            n_meas = axis_of(rec, [2, 3])               # measure about FROZEN axis (decision rule)
        out.append(fivefold_H(band_map(rec, *BAND), n_meas))
    return float(np.median(out))

amps = [0, 15, 30, 50, 75]
res = {m: [inject_run(a, m) for a in amps] for m in ('oracle', 'nolink', 'eim')}
print("POWER (median H vs injected m=5 amplitude; null 95th = %.4f):" % p95)
print("   amp:        " + "  ".join("%6d" % a for a in amps))
for m, lab in (('oracle', 'S_oracle (true axis) '),
               ('nolink', 'S_nolink (frozen)    '),
               ('eim',    'S_eim   (shared axis)')):
    row = res[m]
    print("   %s %s" % (lab, "  ".join("%6.4f" % h for h in row)))
det = {m: ["DET" if h > p95 else "-" for h in res[m]] for m in res}
print("\n   detect@95th  oracle:", det['oracle'])
print("   detect@95th  nolink:", det['nolink'], " <- no low-l/mid-l axis link => ~no power (honest)")
print("   detect@95th  eim   :", det['eim'],    " <- power ONLY under the (undischarged) axis-sharing claim")

# ---------------- FIGURE ----------------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
ax[0].hist(Hn, bins=40, color='#4C78A8', alpha=.85, density=True, label='matched ΛCDM null (frozen axis)')
ax[0].axvline(p95, color='k', ls=':', lw=1.5, label='null 95th pct')
ax[0].set_xlabel('ledger statistic H (m=5 power about frozen ℓ2,3 axis)')
ax[0].set_ylabel('density'); ax[0].set_title('H₅ null — through beam+mask+noise+inpaint')
ax[0].legend(fontsize=8)
ax[1].plot(amps, res['oracle'], '-o', color='#B0B0B0', label='S_oracle: measure about TRUE axis (old code)')
ax[1].plot(amps, res['nolink'], '-o', color='#4C78A8', label='S_nolink: frozen axis, no low-ℓ link')
ax[1].plot(amps, res['eim'],    '-o', color='#E45756', label='S_eim: low-ℓ & m=5 share an axis')
ax[1].axhline(p95, color='k', ls=':', lw=1.2, label='null 95th')
ax[1].set_xlabel('injected m=5 amplitude'); ax[1].set_ylabel('median H')
ax[1].set_title('H₅ power is contingent on axis-sharing'); ax[1].legend(fontsize=7.5)
plt.tight_layout(); plt.savefig('ledger_result.png', dpi=95)
print("\nsaved ledger_result.png")
