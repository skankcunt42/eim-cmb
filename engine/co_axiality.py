"""
co_axiality.py  --  alpha = arccos|n_even . n_odd| co-axiality (l=2 vs l=3), HARDENED.

WHAT CHANGED vs the original:
  (1) MATCHED NULL. alpha is now measured on axes recovered THROUGH the pipeline (beam + mask +
      noise + inpaint), not on ideal full-sky maps. If the mask breaks isotropy, the null median
      moves off 60 deg and the observed-10deg percentile changes -- so the percentile is quoted
      against the matched null, not the analytic 1-cos(alpha).
  (2) FRAME AUDIT BUILT IN (the load-bearing control, CMB-08). We test whether the pipeline itself
      manufactures a preferred axis by asking where the recovered even/odd axes point relative to
      the mask symmetry axis (the 'Galactic pole', z here). A mask-induced axis clustering would
      make ANY co-axiality partly an artifact. This is the control that decides whether an observed
      small alpha is cosmic or frame-driven.

The axis estimator self-test (pure (2,2)/(3,3) -> z) is retained unchanged as an instrument check.
"""
import warnings; warnings.simplefilter('ignore')
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import healpy as hp
import eim_pipeline as E
from eim_pipeline import (LMAX, get_cl, draw_true_alm, observe, make_mask,
                          axis_of, ang, lset_map)

rng = np.random.default_rng(7)
cl = get_cl(LMAX)
mask = make_mask(gal_cut_deg=20.0, apod_deg=5.0)
POLE = np.array([0., 0., 1.])          # mask symmetry axis in these coords ('Galactic pole' analog)
FWHM, NOISE = 5.0, 1.0
NMC = 600

def idx(L, m): return hp.Alm.getidx(LMAX, L, m)

# ================= SELF-TEST (instrument): pure (2,2)/(3,3) planar -> z axis =================
a0 = np.zeros(hp.Alm.getsize(LMAX), complex); a0[idx(2, 2)] = 1.0
e1 = ang(axis_of(a0, [2]), POLE)
a0[:] = 0; a0[idx(3, 3)] = 1.0
e2 = ang(axis_of(a0, [3]), POLE)
print("SELF-TEST: (2,2)->z err=%.1f deg  (3,3)->z err=%.1f deg  [pass if <5]" % (e1, e2))
assert e1 < 5 and e2 < 5, "axis estimator failed self-test"
print("SELF-TEST PASSED\n")

# ================= MATCHED NULL + FRAME AUDIT =================
alpha23 = np.zeros(NMC)
ang_even_pole = np.zeros(NMC)          # frame audit: angle of recovered even axis to mask pole
for i in range(NMC):
    a = draw_true_alm(cl, rng=rng)
    _, rec = observe(a, fwhm_deg=FWHM, mask=mask, noise_uK=NOISE, rng=rng)
    ne = axis_of(rec, [2]); no = axis_of(rec, [3])
    alpha23[i] = ang(ne, no)
    ang_even_pole[i] = ang(ne, POLE)
med = np.median(alpha23)
frac10 = (alpha23 < 10).mean()
print("MATCHED NULL alpha(l2,l3) through beam+mask+noise+inpaint (NMC=%d):" % NMC)
print("   median=%.1f deg   (ideal isotropic 60 deg; deviation = mask/inpaint coupling)" % med)
print("   frac(alpha<10deg)=%.4f   (analytic isotropic 1-cos10 = %.4f)" % (frac10, 1 - np.cos(np.radians(10))))

# frame audit (TWO-SIDED): an isotropic pipeline gives axis-to-pole angle uniform in cos,
# i.e. median 60 deg and CDF F(a)=1-cos a. Deviation EITHER way (pole-clustering OR pole-avoidance)
# is a mask-induced axis anisotropy that any real-axis claim must be compared against.
med_pole = np.median(ang_even_pole)
cospole = np.cos(np.radians(ang_even_pole))                 # ang=arccos|dot| in [0,90] => cos in [0,1]
ks = float(np.max(np.abs(np.sort(cospole) - np.linspace(0, 1, len(cospole)))))   # vs U[0,1] (isotropic)
f_near_pole = (ang_even_pole < 20).mean(); iso20 = 1 - np.cos(np.radians(20))
anisotropic = abs(med_pole - 60) > 6 or f_near_pole > 2 * iso20 or f_near_pole < 0.4 * iso20
print("   FRAME AUDIT (two-sided): median angle(even axis, mask pole)=%.1f deg (iso 60); KS(cos vs U)=%.3f." % (med_pole, ks))
print("                frac within 20deg of pole=%.3f (iso %.3f) -> %s" %
      (f_near_pole, iso20,
       "MASK INDUCES AXIS ANISOTROPY (here: pole-avoidance) -> co-axiality is partly a frame effect"
       if anisotropic else "axis directions consistent with isotropy"))

# ================= observed anchor, quoted against the MATCHED null =================
obs_alpha = 10.0
pct_matched = float((alpha23 <= obs_alpha).mean())
print("\nANCHOR: literature quad-octopole align ~%.0f deg." % obs_alpha)
print("   percentile vs MATCHED null = %.4f   (vs analytic isotropic %.4f)" % (pct_matched, 1 - np.cos(np.radians(obs_alpha))))
print("   NOTE: the known ~10deg alignment is ecliptic/dipole-correlated; the frame audit above,")
print("         not this percentile, is what decides cosmic vs systematic. Real Planck axes plug in here.")

# ================= POWER: impose shared planar z-axis on l=2,3, through the pipeline =================
def _amp(L): return np.sqrt((2 * L + 1) * cl[L])
strengths = [0.0, 0.25, 0.5, 0.75, 1.0]; res_inj = []
for w in strengths:
    got = []
    for _ in range(30):
        a = draw_true_alm(cl, rng=rng)
        for L in (2, 3):
            for m in range(L + 1):
                a[idx(L, m)] *= (1 - w)
            a[idx(L, L)] += w * _amp(L)
        _, rec = observe(a, fwhm_deg=FWHM, mask=mask, noise_uK=NOISE, rng=rng)
        got.append(ang(axis_of(rec, [2]), axis_of(rec, [3])))
    res_inj.append(np.median(got))
    print("   aligned fraction w=%.2f -> median alpha(l2,l3)=%5.1f deg" % (w, np.median(got)))

# ================= FIGURE =================
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
ax[0].hist(alpha23, bins=30, color='#4C78A8', alpha=.85, density=True, label='matched ΛCDM null α(ℓ2,ℓ3)')
aa = np.linspace(0.5, 89.5, 200)
ax[0].plot(aa, np.sin(np.radians(aa)) * (np.pi / 180), 'k--', lw=1.5, label='ideal isotropic  sin α')
ax[0].axvline(10, color='#E45756', lw=2, label='observed ~10° (matched pct=%.3f)' % pct_matched)
ax[0].set_xlabel('α between even & odd axis (deg)'); ax[0].set_ylabel('density')
ax[0].set_title('Co-axiality null (matched) vs known anomaly'); ax[0].legend(fontsize=8)
ax[1].hist(ang_even_pole, bins=30, color='#54A24B', alpha=.8, density=True, label='angle(even axis, mask pole)')
ax[1].plot(aa, np.sin(np.radians(aa)) * (np.pi / 180), 'k--', lw=1.5, label='isotropic sin α')
ax[1].axvline(med_pole, color='#E45756', lw=2, label='median=%.0f°' % med_pole)
ax[1].set_xlabel('angle to mask pole (deg)'); ax[1].set_ylabel('density')
ax[1].set_title('Frame audit: does the mask manufacture an axis?'); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('coaxiality_result.png', dpi=95)
print("\nsaved coaxiality_result.png")
