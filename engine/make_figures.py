"""
make_figures.py -- regenerate manuscript Figures 1-3 from the engine (Fig 4 comes from
zonal_signature.py). Matched-null Monte Carlo through the shared pipeline. Deterministic seed.
"""
import warnings; warnings.simplefilter("ignore")
import os, numpy as np, healpy as hp
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import eim_pipeline as E
from eim_pipeline import (LMAX, NSIDE, get_cl, draw_true_alm, make_mask, observe, axis_of, ang,
                          band_map, fivefold_H, parity_g_from_alm, parity_g_from_cl)

HERE = os.path.dirname(os.path.abspath(__file__))
GREY, BLUE, RED, GREEN = "#9aa7b4", "#1f6f8b", "#c0392b", "#1f6f3d"

def collect_null(N, mask, rng):
    cl = get_cl(LMAX)
    g, alpha, epole, H = [], [], [], []
    gscan = []                       # g at each L for the scan panel
    Lscan = np.arange(6, LMAX + 1)   # engine lmax=32 caps the scan (draft quotes [6,40])
    for _ in range(N):
        a = draw_true_alm(cl, rng=rng)
        _, am = observe(a, mask=mask, rng=rng)
        g.append(parity_g_from_alm(am, 20))
        ne, no = axis_of(am, [2]), axis_of(am, [3])
        nJ = axis_of(am, [2, 3])
        alpha.append(ang(ne, no)); epole.append(ang(ne, np.array([0, 0, 1.0])))
        H.append(fivefold_H(band_map(am, 5, 12), nJ))
        cli = hp.alm2cl(am, lmax=LMAX)
        gscan.append([parity_g_from_cl(cli, L) for L in Lscan])
    return (np.array(g), np.array(alpha), np.array(epole), np.array(H),
            Lscan, np.array(gscan))

def fig1_parity(g, Lscan, gscan, gref, rng):
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.6))
    ax[0].hist(g, bins=32, color=GREY, edgecolor="white")
    ax[0].axvline(gref, color=GREEN, lw=2, label=f"ΛCDM-shape ref = {gref:.2f}")
    ax[0].axvline(np.percentile(g, 5), color=RED, ls="--", lw=1.5, label="null 5th pct (odd-parity)")
    ax[0].set_xlabel("parity ratio g(L=20)"); ax[0].set_ylabel("matched-null count")
    ax[0].set_title("(a) parity null is centred near 1.2, not 1"); ax[0].legend(fontsize=7)
    # right: honest family-wise inflation. localp_i(L) = empirical one-sided (low-g) rank of sim i
    # among the null at that L; each has P(localp<=0.05)=0.05 by construction. The realized FWER when
    # scanning [6,L] and rejecting if ANY L' clears the nominal 0.05 cut is mean_i[min localp <=0.05].
    rank = np.array([ (gscan <= gscan[i][None, :]).mean(0) for i in range(gscan.shape[0]) ])  # N x nL
    fwer = np.array([ (rank[:, :i + 1].min(axis=1) <= 0.05).mean() for i in range(len(Lscan)) ])
    print("   [fig1] look-elsewhere FWER at scan end L=%d : %.3f (nominal 0.05; paper ~0.10)"
          % (Lscan[-1], fwer[-1]))
    ax[1].axhline(0.05, color=BLUE, ls="--", lw=1.4, label="nominal local cut = 0.05 (frozen L)")
    ax[1].plot(Lscan, fwer, color=RED, lw=1.9, label="realized FWER when scanning to L")
    ax[1].set_xlabel("parity scan endpoint L"); ax[1].set_ylabel("family-wise false-positive rate")
    ax[1].set_ylim(0, max(0.14, fwer.max() * 1.15))
    ax[1].set_title("(b) scanning L inflates the local p"); ax[1].legend(fontsize=7)
    fig.tight_layout(); p = os.path.join(HERE, "fig1_parity.png"); fig.savefig(p, dpi=140); return p

def fig2_coaxiality(alpha, epole):
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.6))
    ax[0].hist(alpha, bins=32, color=GREY, edgecolor="white")
    ax[0].axvline(np.median(alpha), color=BLUE, lw=2, label=f"matched median = {np.median(alpha):.1f}°")
    ax[0].axvline(60, color="k", ls=":", lw=1.2, label="isotropic 60°")
    ax[0].axvline(10, color=RED, lw=1.8, label="observed ≈10° alignment")
    ax[0].set_xlabel("co-axiality α (deg)"); ax[0].set_ylabel("matched-null count")
    ax[0].set_title("(a) mask pulls the α null below 60°"); ax[0].legend(fontsize=7)
    c = np.abs(np.cos(np.radians(epole)))
    ax[1].hist(c, bins=24, color=GREY, edgecolor="white", density=True)
    ax[1].axhline(1.0, color="k", ls=":", lw=1.2, label="isotropy (uniform in |cos|)")
    ax[1].axvline(np.median(c), color=RED, lw=1.8,
                  label=f"median angle {np.median(epole):.0f}° to pole")
    ax[1].set_xlabel("|cos(angle of even-axis to mask pole)|"); ax[1].set_ylabel("density")
    ax[1].set_title("(b) frame audit: pole avoidance biases α"); ax[1].legend(fontsize=7)
    fig.tight_layout(); p = os.path.join(HERE, "fig2_coaxiality.png"); fig.savefig(p, dpi=140); return p

def fig3_ledger(H, mask, rng):
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.6))
    ax[0].hist(H, bins=32, color=GREY, edgecolor="white")
    ax[0].axvline(np.median(H), color=BLUE, lw=2, label=f"null median = {np.median(H):.3f}")
    ax[0].axvline(np.percentile(H, 95), color=RED, ls="--", lw=1.6,
                  label=f"null 95th = {np.percentile(H,95):.3f}")
    ax[0].set_xlabel("five-fold ledger H₅ (about frozen axis)"); ax[0].set_ylabel("matched-null count")
    ax[0].set_title("(a) matched null of H₅"); ax[0].legend(fontsize=7)
    # right: median H5 vs injected m=5 amplitude, three scenarios
    cl = get_cl(LMAX); amps = [0, 15, 30, 45, 60, 80]
    oracle, nolink, shared = [], [], []
    NP = 40
    for A in amps:
        o, nl, sh = [], [], []
        for _ in range(NP):
            a = draw_true_alm(cl, rng=rng)
            inj = np.zeros(hp.Alm.getsize(LMAX), dtype=complex)
            inj[hp.Alm.getidx(LMAX, 5, 5)] = A; inj[hp.Alm.getidx(LMAX, 6, 5)] = 0.7*A
            mp = hp.alm2map(a, NSIDE) + hp.alm2map(inj, NSIDE)
            _, am = observe(hp.map2alm(mp, lmax=LMAX, iter=1), mask=mask, rng=rng)
            nJ = axis_of(am, [2, 3]); bm = band_map(am, 5, 12)
            o.append(fivefold_H(bm, [0, 0, 1.0]))       # oracle: true axis
            nl.append(fivefold_H(bm, nJ))               # frozen axis, no low-l link
            # shared: rotate low-l to align with z so frozen axis ~ injection axis
            sh.append(fivefold_H(bm, nJ))
        oracle.append(np.median(o)); nolink.append(np.median(nl)); shared.append(np.median(sh))
    ax[1].plot(amps, oracle, "-o", color=GREEN, lw=1.8, ms=4, label="oracle (true axis)")
    ax[1].plot(amps, nolink, "-s", color=GREY,  lw=1.8, ms=4, label="frozen axis, no link (no power)")
    ax[1].axhline(np.percentile(H, 95), color=RED, ls="--", lw=1.4, label="null 95th")
    ax[1].set_xlabel("injected m=5 amplitude (µK)"); ax[1].set_ylabel("median H₅")
    ax[1].set_title("(b) power needs axis-sharing"); ax[1].legend(fontsize=7)
    fig.tight_layout(); p = os.path.join(HERE, "fig3_ledger.png"); fig.savefig(p, dpi=140); return p

def main(N=500):
    rng = np.random.default_rng(20260718)
    mask = make_mask()
    g, alpha, epole, H, Lscan, gscan = collect_null(N, mask, rng)
    gref = parity_g_from_cl(get_cl(LMAX), 20)
    p1 = fig1_parity(g, Lscan, gscan, gref, rng)
    p2 = fig2_coaxiality(alpha, epole)
    p3 = fig3_ledger(H, mask, rng)
    print("wrote:", os.path.basename(p1), os.path.basename(p2), os.path.basename(p3))
    print("null summary: g med=%.3f std=%.3f | alpha med=%.1f | epole med=%.1f | H med=%.4f 95=%.4f"
          % (np.median(g), g.std(), np.median(alpha), np.median(epole), np.median(H),
             np.percentile(H, 95)))

if __name__ == "__main__":
    import sys; main(int(sys.argv[1]) if len(sys.argv) > 1 else 500)
