# EIM–CMB Parity Engine — Predeclaration & Rulings

**Filed 2026-07-18 · before any run on real Planck data · radical-doubt frame**

This document consolidates the closure trajectory for the EIM/CMB refractive-index
hypothesis into a predeclared, four-statistic test engine plus the standing rulings.
Everything here is filed **before** touching real data. The predeclared choices
(axis source, ℓ-ranges, statistics, decision rules) may **not** be changed post-hoc.
A predeclaration that is edited after seeing the data is void.

---

## 0. One-line status

EIM survives first contact **only** as a finite, predeclared CMB parity / co-axiality /
five-fold-ledger hypothesis — **not** as a derivation of General Relativity.

---

## 1. Standing rulings (do not silently undo)

- **No β(n̂) polarization reconstruction** until EIM predeclares a **sign** or an **axis**
  (see `EIM_CMB_Polarization_Predeclaration_Memo.docx`). An unconstrained reconstruction is a
  confirmation-bias net: over any sky it returns a strongest axis and a nonzero amplitude.
- **The dodecahedral / S² backbone is L4-conditional.** The construction D→X→Γ is open;
  the flag validator is **not** evidence EIM produced the input; A11 retracted (D ≠ Γ).
- **Observer / quotient are frame/readout operations, not primitive ontology.**
- **GR derivation is undefined-within-processor.** ℤ₂ ⇏ SO(3,1); anchor-count ⇏ T_ab;
  G/c⁴ is not produced. Do **not** cite G = 8πT as derived. Curvature = anholonomy is a real
  kernel but is gated on an explicit ℤ₂ → SO(3,1) coarse-graining that does not yet exist.

---

## 2. The predeclared engine  (temperature scalar  n(Ω) − 1 = ΔT(Ω)/T₀)

**Frozen-axis protocol.** The reference axis **n̂_J is fixed from ℓ ∈ {2,3} ONLY** (the co-axis),
before any other statistic is computed. No statistic may re-fit the axis. This is what makes the
ledger test non-post-hoc.

**Statistics** — all compared to **matched ΛCDM Monte Carlo** with the same beam, mask, noise,
smoothing, and inpainting; **never** to ideal full-sky isotropy:

| Symbol | Definition | ℓ-range | Role |
|---|---|---|---|
| g_L | P⁺/P⁻ even/odd multipole power ratio | ℓ ≤ L, **L predeclared** | parity (do NOT scan L to minimise p) |
| α | arccos \|n̂_even·n̂_odd\| co-axiality | ℓ=2 vs ℓ=3 | known-anomaly baseline |
| A_antipodal | corr( n(Ω), n(−Ω) ) | low-ℓ | ≈ g_L by even/odd↔antipodal identity |
| **H₅** | fractional **m=5** azimuthal power about **n̂_J** | ℓ ∈ [5,12] | **primary EIM-specific differential test** |

**Decision rule.** EIM gains a real hit **iff** the joint { g_L, α, H₅ } beats matched sims across
≥ 2 component-separation maps (SMICA / NILC / SEVEM / Commander) and ≥ 2 masks, look-elsewhere
corrected, **and** H₅ locks to the **same** n̂_J across independent pipelines. Absent that, close
as known-anomaly / interpretation.

---

## 3. Prior weights (honest, before the run)

- **g_L** — known ~2–3σ odd-parity preference, but look-elsewhere-marginal (z slides with L).
- **α** — known ~1.5% quad-octopole co-axiality, but ecliptic/dipole-contaminated.
- **A_antipodal** — ≈ g_L; report but do not double-count.
- **H₅** — **the real test.** Not a catalogued anomaly, not ecliptic-degenerate by construction,
  and **provably independent of α** (measured corr(H₅, α) = −0.03 under the null). Prior to note:
  the Poincaré dodecahedral *topology* was searched by Planck matched-circles and **disfavoured** —
  H₅ is a distinct **parity/closure-cell** signature (five-fold power locked to the frozen ℓ=2,3
  axis), not a global-topology claim. Nobody has run it.

---

## 4. Engine validation status (this bundle)

No real Planck data was accessible in the build environment; every pipeline below is
**self-tested and null-validated on matched ΛCDM simulations**. Plug in real low-ℓ a_ℓm to
obtain real-sky numbers — no rebuild required.

- `parity_test.py` — g_L = P⁺/P⁻; CAMB ΛCDM theory Cℓ; cosmic-variance MC null; power via
  even-ℓ suppression. Shows the z-vs-ℓ_max look-elsewhere slide explicitly.
- `co_axiality.py` — convention-free **power-tensor axis estimator** (self-test: pure (2,2)/(3,3)
  → z, 0.0° error); null α(ℓ2,ℓ3) matches isotropic (median 60°); power via imposed aligned fraction.
- `check_cdf.py` — high-N (8000) CDF vs analytic 1−cos α; **confirms no floor/selection**
  (max|CDF−(1−cos α)| = 0.010); exact MC percentile of the 10° co-axiality = 0.0154.
- `build_birefringence.py` — isotropic-β EB estimator; self-test recovers injected β (0.30→0.301,
  unbiased); Fisher/MC sensitivity σ_β; detectability vs noise. Flags the α–β calibration degeneracy.
- `build_ledger.py` — **H₅** five-fold (m=5) power about the frozen axis; self-test (5-fold about z
  detected 6:1); ΛCDM null; **independence corr(H₅,α) = −0.03**; power via injected m=5.

---

## 5. Final status ledger

- **Exact graph fact:** τ²=1 antipodal closure; even/odd (ℤ₂) split; 6 face-pair decomposition;
  the five-fold azimuthal structure about a chosen axis.
- **Processor-specific postulate:** an unresolved single Janus circuit favours a stable
  odd / parity / co-axial / five-fold residue.
- **Interpretive mapping:** CMB scalar anisotropy as n_cosmic(Ω) − 1.
- **Speculative cosmological claim:** the observed CMB anomalies *are* that residue.
- **Undefined-within-processor:** G/c⁴; T_ab from scalar counts; SO(3,1) holonomy from ℤ₂ flips;
  the dodecahedral / S² backbone from EIM primitives.

---

## 6. The one experiment that matters

Freeze n̂_J from ℓ=2,3 only → compute **H₅** about it in ℓ∈[5,12] → require it to beat matched
ΛCDM across SMICA/NILC/SEVEM/Commander × two masks → require the **same axis** to carry the excess
in independent pipelines. If H₅ locks and beats the matched null, EIM has a real, non-post-hoc,
Lorentz-safe hit that is not a known anomaly. If not, the arc closes as interpretation. **Do not
drag Einstein's equation into the proof until H₅ survives a matched null.**

---

## 7. Manifest

```
PREDECLARATION_AND_RULINGS.md            this file
EIM_CMB_Polarization_Predeclaration_Memo.docx   sign/axis ruling (β route → interpretive)
parity_test.py                           g_L parity statistic + null + power
co_axiality.py                           power-tensor axis estimator + null + power
check_cdf.py                             co-axiality CDF check (no-floor proof)
build_birefringence.py                   isotropic-β EB estimator + null + sensitivity
build_ledger.py                          H₅ five-fold ledger + null + independence + power
```

Missing inputs to close the trajectory: real Planck low-ℓ a_ℓm (or component-separated maps +
masks), and a signature on this predeclaration before the run.
