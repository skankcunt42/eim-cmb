# Predeclaration addendum — the zonal five-fold signature (§5 of the closure paper)

**Dated:** 2026-07-18 · **Status:** frozen before real-data contact with the zonal statistic.
**Supersedes for the framework test:** the five-fold ledger H₅, which the closure paper (§4)
reclassifies as an *anomaly search*, not a test of the dodecahedral hypothesis. This addendum
freezes the observable the hypothesis actually predicts (paper §5) so that a real-data run is
look-elsewhere–protected in the same spirit as the original predeclaration.

This document fixes definitions and the decision rule. It does not report any data value. The
implementing code is `zonal_signature.py` (statistic + matched null) and `run_closure.py` (the
≥2-map × ≥2-mask driver); the group-theoretic basis is verified in `selection_rules.py`.

## 1. Motivation (fixed, not re-derived here)
By Theorem 2 (k-fold azimuthal selection, verified to machine precision in `selection_rules.py`),
a field invariant under rotation by 2π/5 about an axis n̂ has azimuthal content only in m ≡ 0
(mod 5); hence about a dodecahedral five-fold axis the multipoles ℓ = 2,3,4 are **zonal** (m = 0)
and the first non-axial power is at m = 5, ℓ = 5. A genuine five-fold axis is therefore *zonal*,
whereas the observed axis-of-evil morphology is *sectoral* (planar octopole, m = ±3). The zonal
signature is the falsifiable discriminant.

## 2. Statistic definitions (frozen)
Let a map's recovered harmonic coefficients (after the matched pipeline, §4) be a_ℓm.

- **Zonal axis** n̂_Z ≔ the **largest-eigenvalue** eigenvector of the power tensor
  M_ij = Σ_p T(p)² v_i(p) v_j(p), where T is the map synthesized from ℓ ∈ {2,3,4} and v(p) are pixel
  unit vectors. (The sectoral estimator used for the axis of evil is the *smallest*-eigenvalue
  eigenvector; on a zonal pattern the two are perpendicular — verified in `selection_rules.py`.)
  n̂_Z is frozen from ℓ ∈ {2,3,4}; **no statistic below is permitted to re-fit it.**
- **zonal_low** ≔ mean over ℓ ∈ {2,3,4} of the m = 0 power fraction about n̂_Z, i.e.
  |a_ℓ0|² / Σ_{m=−ℓ}^{ℓ} |a_ℓm|² evaluated in the frame with n̂_Z → ẑ.
- **m5_mid** ≔ mean over ℓ ∈ {5,6} of the m = 5 power fraction about n̂_Z, i.e.
  2|a_ℓ5|² / Σ_m |a_ℓm|² in the same frame.

The ℓ-sets {2,3,4} (axis + zonal_low) and {5,6} (m5_mid) are **fixed**; no scan over ℓ-ranges is
permitted. Frame rotation uses the verified Euler convention in `zonal_signature.py`.

Geometry check (dimensionless, from `zonal_signature.py`): a genuine zonal five-fold sky →
(zonal_low, m5_mid) ≈ (1, 1) about any axis; a realistic sectoral axis-of-evil sky and the ΛCDM
null both sit near (0.26, 0.17). The genuine signature occupies a region of the plane the null
leaves empty; the observed morphology does not.

## 3. Null and thresholds (frozen)
The decision null is the **matched pipeline null** — ΛCDM skies carried through the identical
beam (5°) + apodized symmetric mask + homogeneous white noise + iterative-harmonic inpainting +
harmonic recovery used on the data. The per-(map, mask) decision thresholds are the **95th
percentiles** of the matched-null zonal_low and m5_mid.

> **Do not** use the full-sky isotropic null for the decision. The paper's Figure-4 null (0.26,
> 0.17) tracks the *full-sky* null; the matched pipeline suppresses mid-ℓ m = 5 power, moving the
> matched-null m5_mid centre to ≈ 0.12. The framework test is decided against the matched null.

## 4. Decision rule (frozen)
The **zonal five-fold signature is detected** iff **both** zonal_low **and** m5_mid exceed their
matched-null 95th percentile (one-sided, high tail), **jointly**, for **≥ 2 component-separated
maps AND ≥ 2 masks**, with the frozen zonal axes agreeing across those maps to **< 10°**
(common-axis requirement), after look-elsewhere correction. Absent that, the result is reported as
**consistent with ΛCDM / no zonal five-fold signature** — a null that *counts against* the
dodecahedral reading rather than being absorbed as a pre-existing anomaly.

## 5. Maps, masks, corrections (frozen)
- **Maps:** SMICA, NILC, SEVEM, Commander (Planck 2018 PR3 component-separated, Nside = 64, µK,
  monopole+dipole removed) from `fetch_planck_streamlined.ipynb`.
- **Masks:** ≥ 2, both real. Primary = Planck common intensity mask (f_sky≈0.78); secondary =
  Galactic-plane **GAL070** (`HFI_Mask_GalPlane-apo2_2048_R2.00.fits`, field 3, f_sky≈0.70), a
  genuinely different construction and sky fraction. Both are produced by
  `fetch_planck_streamlined.ipynb`; `run_closure.py` falls back to a synthetic second cut only if
  the real GAL070 mask is absent, and prints a warning when it does.
- **Look-elsewhere:** no ℓ-range scan (§2). The zonal axis is *estimated* by the fixed estimator,
  not searched. Any auxiliary L-scan (parity) is reported only to expose inflation, never as the
  decision value.
- **Inpainting:** iterative-harmonic here (applied identically to null and data); swap for the
  survey map-maker's inpainter at submission.

## 6. Epistemic label (carried from the paper)
A positive zonal signature would support the processor-specific postulate that a single closure
cell selects one five-fold face-pair (breaking I_h to a five-fold subgroup); it would **not**
derive that postulate from first principles, and the dodecahedral backbone remains a conditional
realization. H₅ rows, if reported, are anomaly documentation, not a framework test.
