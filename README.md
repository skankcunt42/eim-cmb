# EIM–CMB

Reproducible engine for a preregistered, matched-null test of the dodecahedral / five-fold
(EIM closure) hypothesis for the large-angle CMB temperature anomalies — plus the manuscript that
reports it. Everything runs **in the browser**; nothing needs to be installed locally.

## Run it in the browser (no local Python)

1. Open **`colab/RUN_EVERYTHING.ipynb`** in Google Colab (badge below, or File → Open notebook →
   GitHub → this repo).
2. Runtime → **Run all**.

That clones this repo, installs the dependencies, runs the data-independent self-tests
(Theorems 1–2, estimator checks), regenerates the paper's §III numbers, rebuilds all four figures,
**downloads the real Planck 2018 maps** (Colab has open internet — the fetch that a locked sandbox
can't do runs here), and executes the full closure driver over four maps × two masks.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GITHUB_USER/eim-cmb/blob/main/colab/RUN_EVERYTHING.ipynb)

> Replace `GITHUB_USER` in the badge and in the first Colab cell with the account this repo lives under.

## What's here

```
engine/     the computational core (import-and-run Python modules)
  eim_pipeline.py          shared beam/mask/noise/inpaint matched pipeline
  parity_test.py           g_L parity ratio, matched null, look-elsewhere
  co_axiality.py           alpha co-axiality + frame audit
  build_ledger.py          H5 five-fold ledger, three-scenario power
  check_cdf.py             full-sky no-floor instrument check
  selection_rules.py       Theorems 1 & 2 + estimator self-tests (machine precision)
  zonal_signature.py       the section-5 zonal observable + Figure 4
  reproduce_paper_numbers.py   independent regeneration of the section-3 numbers
  make_figures.py          Figures 1–3
  run_closure.py           full-closure driver: >=2 maps x >=2 masks
  stream_fetch.py          streaming Planck fetch (download -> extract T -> downgrade -> compress)
notebooks/  fetch_planck_streamlined.ipynb   step-by-step data acquisition
colab/      RUN_EVERYTHING.ipynb             one-click end-to-end run
docs/       predeclaration, hardening log, testing/run instructions
paper/      EIM_CMB_closure_paper_FINAL.docx  the manuscript (PRD full article)
figures/    the four manuscript figures (PNG)
```

## How we work on this

- **GitHub (this repo)** is the permanent, versioned home — the single source of truth.
- **Colab** is the compute surface — browser only, real internet, zero install.
- Results (`*_nside64_uK.fits`, `*.npz`, figures) are *regenerated*, never committed (see
  `.gitignore`); commit code and docs, not data.

## Status

The manuscript (`paper/`) is a preregistration + methodology + group-theoretic obstruction paper
whose conclusions do not depend on a new map analysis. Section VI (the real-data run) is a scaffold
with placeholder tables; running `colab/RUN_EVERYTHING.ipynb` on the real maps populates it. The
group-theoretic obstruction (Theorems 1–2) is verified to machine precision and is data-independent.

See `docs/PREDECLARATION_ZONAL_ADDENDUM.md` for the frozen §V statistic and `docs/TESTING_README.md`
for the module-by-module validation.
