# EIM–CMB Parity Engine — Handoff Package

Hardened, predeclared, matched-null test engine for the EIM/CMB refractive-index hypothesis.
Everything needed to run the frozen-axis g_L / α / H₅ protocol is here.

## What's in this package
- `PREDECLARATION_AND_RULINGS.md` — the frozen protocol (axis source, ℓ-ranges, statistics,
  decision rule) and standing rulings. Read first.
- `HARDENING_LOG.md` — what was hardened and the substantive findings (axis–power fix; parity
  g=1-reference flaw → true ΛCDM-shape reference ≈1.13; co-axiality mask bias 60°→54°; H₅ power is
  contingent on the undischarged axis-sharing bridge). Read second.
- `RUN_INSTRUCTIONS.md` — exact steps to run on real Planck data.
- `eim_pipeline.py` — shared matched pipeline (beam+mask+noise+inpaint). Imported by the rest.
- `parity_test.py`, `co_axiality.py`, `build_ledger.py`, `check_cdf.py` — the four self-testing,
  matched-null statistic engines (run on ΛCDM sims out of the box, no data needed).
- `fetch_planck.py` — downloads official Planck 2018 maps + common mask, downgrades to Nside=64,
  writes ~1 MB files. Run where the Planck archive is reachable (your machine).
- `run_real.py` — executes the frozen-axis protocol on the real maps produced by fetch_planck.py.

## Quick start
```
pip install healpy camb astropy numpy scipy matplotlib
# engine self-tests + matched-null demos (no data required):
python parity_test.py       # parity g_L, matched null, look-elsewhere
python co_axiality.py       # co-axiality α, matched null, frame audit
python build_ledger.py      # H₅ five-fold ledger, three-scenario power
python check_cdf.py         # full-sky no-floor instrument check
# real data (needs Planck-archive network access):
python fetch_planck.py      # -> smica/nilc_nside64_uK.fits, commonmask_nside64.fits
python run_real.py          # -> real-sky g_L, α, H₅ with matched-null p-values
```

## One caveat that matters for interpretation
H₅ is a test *of EIM* only under the axis-sharing assumption the framework does not yet derive
(HARDENING_LOG finding #1). Absent that bridge, a surviving H₅ is an anomaly loosely motivated by
EIM, not a confirmation of it. Every hardening correction moved significance down or revealed a
frame dependence — consistent with the predeclaration's honest prior.
