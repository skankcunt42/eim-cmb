# EIM–CMB Parity Engine — Hardening Log

**Filed 2026-07-18 · post-predeclaration engine hardening (still before any real-data run).**

Scope of this pass (chosen target: *harden the engine* so the predeclared test is
publication-grade before any real Planck run). Nothing here touches the predeclared
*choices* (axis source, ℓ-ranges, statistics, decision rule) — those remain frozen. What
changed is that every statistic is now computed through a **matched simulation pipeline**
and every claimed detectability is now measured with the **decision-rule statistic**, not an
optimistic proxy. Two of the changes are not tidying: they move the numbers.

---

## 0. New shared module

`eim_pipeline.py` — one pipeline used identically for null and (future) data:
CAMB ΛCDM Cℓ → beam (5° FWHM) → map → homogeneous white noise → **mask** (Galactic cut,
apodized) → **iterative harmonic inpainting** → recovered aℓm. The null therefore *contains*
the mask/beam/inpaint coupling instead of ignoring it. `observe()` returns recovered aℓm and
is a drop-in for real data (swap `make_mask()` for the real analysis mask; swap `inpaint()`
for the map-maker's inpainter).

## 1. Defect fixed — axis–power mismatch in the H₅ ledger (`build_ledger.py`)

The original power test injected m=5 about *z*, froze `nfroz = axis_of(alm,[2,3])`, then measured
H about *z* — i.e. about the **true** signal axis, not the **frozen** axis the decision rule uses.
That silently overstated detectability. Fixed: H is measured about the frozen axis. The old
number is kept only as a labelled **oracle** upper bound. Three scenarios now make the
consequence explicit (matched null 95th percentile = 0.078):

| injected m=5 amp | 0 | 15 | 30 | 50 | 75 |
|---|---|---|---|---|---|
| S_oracle (measure about TRUE axis — old code) | .033 | .041 | .060 | **.099** | **.124** |
| S_nolink (frozen axis, low-ℓ random) | .047 | .047 | .044 | .049 | .045 |
| S_eim (low-ℓ **and** m=5 share an axis) | .030 | .034 | .055 | .076 | **.104** |

**Reading.** S_nolink never clears the null: with no low-ℓ↔mid-ℓ axis link, the frozen-axis H₅
has *no power*. Power returns only in S_eim, where the ℓ=2,3 axis is strongly co-axial
(aligned fraction 0.9 → frozen axis ~5° from truth, i.e. the ~10°-co-axiality anomaly level)
*and* the five-fold genuinely shares it — and even then it sits below the oracle (frozen-axis
penalty). **This is the axis-sharing bridge, quantified:** H₅ is a test *of EIM* only under an
assumption EIM does not yet derive (polarization memo: τ is structurally axis-free; cellulation
orientation is gauge). Absent that bridge, H₅ is an EIM-*motivated* anomaly hunt, not a test of
the framework. A co-axiality of ~13° (aligned fraction 0.75) already erodes it to borderline.

## 2. Defect fixed — parity reference and null (`parity_test.py`)

Two corrections, both reducing any apparent odd-parity anomaly:

- **Wrong null center.** The ΛCDM-*shape* expectation of g = Σ_even Dℓ / Σ_odd Dℓ (2≤ℓ≤20) is
  **≈1.13**, not 1 — the even-ℓ sum carries the large ℓ=2 Sachs-Wolfe term. The original engine's
  "g=1 expectation" was the flaw; measuring a deficit from 1.0 inflates significance.
- **Wrong variance.** Full-sky σ(g₂₀)≈0.198 → **masked σ≈0.27** (≈36% wider). Matched-null mean
  ≈1.17–1.20.
- **Look-elsewhere quantified.** Scanning L∈[6,40] to minimize p inflates a nominal local
  p=0.05 to a global **p≈0.10 (×2.0)**. Decision number = local p at the frozen L₀=20; the scan
  value is shown only to expose the inflation.

## 3. Defect fixed — co-axiality null + frame audit (`co_axiality.py`)

- Matched-null median α(ℓ2,ℓ3) = **54.1°**, not the isotropic 60°: the mask biases α *low*. So the
  observed ~10° sits at **p=0.020 against the matched null**, vs 0.0152 analytic — *less*
  significant, not more.
- **Frame audit (built in, two-sided).** Recovered even/odd axes are not isotropic under the
  pipeline: median axis-to-pole angle 73.6° (iso 60°), near-pole fraction 0.007 vs isotropic
  0.060, KS(cos vs U[0,1]) = 0.26. The mask *avoids the pole* and confines axes to the
  equatorial band — which is the mechanism that biases α low. **Co-axiality is therefore partly a
  frame effect**; the frame audit (CMB-08), not the raw percentile, is what would decide cosmic
  vs systematic on real data.
- `check_cdf.py` retained unchanged as a *full-sky instrument check* (median 59.93°, max|CDF −
  (1−cos α)| = 0.006, no floor). It now serves as the isotropic control against which the masked
  54.1° is read.

## 4. Verification

- Axis estimator self-test: (2,2)→z, (3,3)→z at 0.0° error. ✓
- Five-fold estimator self-test: H(z)/H(x) = 5.5 (median of 20 injected realizations). ✓
- **Null calibration:** null-vs-null H₅ p-values are uniform (KS D=0.070, p=0.10; mean p=0.53;
  false-positive rate 0.040 vs nominal 0.05). The matched null does not self-generate detections. ✓
- All four figures regenerate from the hardened scripts.

## 5. Residual limits (explicit — these gate a real-data submission)

1. **Synthetic mask.** A symmetric apodized Galactic cut, not the Planck common mask. Real mask
   → replace `make_mask()`; frame-audit anisotropy will differ and must be re-measured.
2. **Simple inpainting.** Iterative harmonic ("cooling"), applied identically to null and data so
   its bias is absorbed — but it is not a constrained Gaussian realization. For submission, use
   the map-maker's inpainter.
3. **No real data.** All numbers are matched-null self-tests. Plug real low-ℓ aℓm into `observe()`.
4. **H₅-as-EIM-test is still gated on the axis-sharing bridge** — a *theory* task (finding #1),
   not a code task. Until that bridge is discharged or H₅ is formally reclassified as
   EIM-motivated anomaly-hunt, a surviving H₅ is an anomaly, not a confirmation of EIM.

## 6. Net effect on the closure

Every hardening correction moved significance **down** or revealed a **frame dependence**: the
parity anomaly loses its g=1 reference and gains mask variance; the co-axiality loses ~6° of null
median to the mask and is shown to be partly frame-driven; H₅'s only route to power runs through
an undischarged bridge. Consistent with the predeclaration's honest prior. The engine is now
matched-null, decision-rule-consistent, look-elsewhere-aware, and calibrated — i.e. ready to
return a *trustworthy* number when real Planck low-ℓ aℓm are supplied, whichever way it falls.
