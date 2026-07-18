# EIM–CMB testing layer (built on the closure-paper baseline)

New modules that turn the closure paper into a runnable, self-validating test suite. They sit on
top of the existing hardened engine (`eim_pipeline.py` + the four statistic modules) and were each
validated against the paper's own stated numbers.

## Modules

- **`selection_rules.py`** — verifies the group-theoretic backbone (§4) and Appendix-A self-tests.
  Theorem 1 (icosahedral selection): P(ℓ)/P_max ≤ 1e-27 for ℓ=1–5, O(1) at ℓ=6,10,12 (direct-aₗₘ,
  cleaner than the paper's 1e-13 pixelized claim). Theorem 2 (k-fold selection): C₅/C₃/C₂ leak
  ≤ 1e-22 into forbidden m. Sectoral↔zonal estimators are 0°/90° on Y₂₂/Y₃₃ vs Y₂₀/Y₃₀.
  *Note:* the five-fold z-to-x ratio is amplitude-dependent (25× pure injection; ≈6 at 80 µK on a
  ΛCDM background); the paper's 5.5 is an injection-on-background number — state the amplitude in
  Appendix A for reproducibility.

- **`zonal_signature.py`** — the §5 observable the hypothesis actually predicts: zonal axis
  estimator + (zonal_low, m5_mid) + matched null. Validated: genuine template → (1,1) about any
  axis; realistic sectoral sky and null both in the low corner; full-sky null (0.260, 0.170) =
  paper (0.26, 0.17). *Finding:* the paper's Fig-4 null is the **full-sky** null; the matched
  pipeline null is (0.273, 0.117) — the decision-relevant one. Writes `zonal_signature_fig4.png`.

- **`reproduce_paper_numbers.py`** — regenerates §3 / Table-1 numbers with MC error bars. At N=400,
  11/11 within tolerance (g ref 1.126; masked std 0.275 vs 0.27; α median 55° vs 54°; frame audit
  74° vs 73.6°; H₅ null 0.046/0.080 vs 0.047/0.078). The H₅ null-vs-null calibration (KS/FPR) needs
  N ≳ 2000 to match the paper's tighter 0.07/0.04.

- **`run_closure.py`** — full-closure driver over ≥2 maps × ≥2 masks, reporting both the legacy
  anomaly set {g, α, H₅} and the §5 framework test {zonal_low, m5_mid}, with per-(map,mask)
  matched-null p-values, frozen-axis cross-map stability, and the two decision rules. Reads the
  four `*_nside64_uK.fits` + common mask if present; otherwise runs a labelled ΛCDM dry-run.

- **`PREDECLARATION_ZONAL_ADDENDUM.md`** — dated freeze of the zonal statistic (definitions,
  matched-null thresholds, ≥2-map × ≥2-mask decision rule, look-elsewhere policy). Read before any
  real-data contact.

## Order of operations
```
python selection_rules.py            # data-independent; verifies Thm 1 & 2 + estimator self-tests
python zonal_signature.py            # §5 statistic self-test + Fig 4 (simulations only)
python reproduce_paper_numbers.py 800  # referee-check §3 numbers (raise N to tighten)
# then, once fetch_planck_streamlined.ipynb has produced the four maps + mask:
python run_closure.py 600            # frozen-axis closure across 4 maps x >=2 masks
```

## Masks (≥2, both real)
The streamlined fetch now produces **two real Planck masks**: `commonmask_nside64.fits` (common
intensity mask, f_sky≈0.78) and `gal070mask_nside64.fits` (Galactic-plane GAL070 from
`HFI_Mask_GalPlane-apo2_2048_R2.00.fits`, field 3, f_sky≈0.70). `run_closure.py` uses both when
present and falls back to synthetic cuts (with a printed warning) only if either is missing. To
vary the second mask, set `MASK2_FIELD=4` in the notebook for GAL080 (f_sky≈0.80).

## Open items for submission
1. Swap the iterative-harmonic inpainter for the survey map-maker's inpainter.
2. Fix Appendix A: state the m=5 injection amplitude behind the 5.5 ratio; and reconcile the
   Fig-4 null (full-sky) with §3's matched-null convention, or label Fig 4 as an isotropic-null
   illustration with the matched-null decision thresholds overlaid.
3. (Optional) Add a third mask tier (e.g. GAL080) for a mask-robustness gradient.
