"""
selection_rules.py -- numerical verification of the group-theoretic backbone of the closure
paper (Theorems 1 and 2) and the estimator self-tests of Appendix A.

Everything here is DATA-INDEPENDENT and (for the selection rules) deterministic: it either
reproduces the paper's claimed machine-precision facts or it does not. No CMB data required.

Theorem 1 (icosahedral selection): a function invariant under the full icosahedral group has
    vanishing multipoles at l=1..5; lowest non-trivial invariant at l=6 (then 10,12). We build a
    sum of unit sources at the 12 dodecahedral face-centre directions (= icosahedron vertices) and
    evaluate the a_lm DIRECTLY from spherical harmonics (no healpix pixelization, so the vanishing
    is exact to floating point, matching the paper's <=1e-13 claim).

Theorem 2 (k-fold azimuthal selection): a function invariant under rotation by 2*pi/k about n has
    azimuthal content only in m == 0 (mod k). We build C_k-symmetric functions about z and show the
    surviving m are exactly the multiples of k, then read off the l-specific corollaries the paper
    uses (5-fold -> l=2,3,4 zonal; 3-fold -> octopole m=+-3 allowed, quadrupole m=0; 2-fold ->
    quadrupole m=+-2 allowed, octopole m=+-3 forbidden).

Appendix A estimator self-tests: sectoral estimator (smallest-eigval power-tensor eigenvector)
    recovers pure Y_{2,2}/Y_{3,3} axes at 0 deg and MISSES zonal axes by 90 deg; zonal estimator
    (largest-eigval eigenvector) does the reverse; the five-fold H5 statistic is far stronger about
    the true m=5 axis than about a perpendicular one.
"""
import warnings; warnings.simplefilter("ignore")
import numpy as np
import healpy as hp

# ---- spherical harmonics (support both modern and legacy scipy signatures) ----
try:
    from scipy.special import sph_harm_y
    def _Ylm(l, m, theta, phi):
        return sph_harm_y(l, m, theta, phi)          # (l, m, theta[polar], phi[azimuth])
except Exception:                                    # legacy scipy
    from scipy.special import sph_harm
    def _Ylm(l, m, theta, phi):
        return sph_harm(m, l, phi, theta)            # (m, l, phi[azimuth], theta[polar])

GOLDEN = (1.0 + np.sqrt(5.0)) / 2.0


# ================================================================ Theorem 1
def icosahedron_vertices():
    """12 icosahedron vertices = 12 dodecahedral face-centre directions (unit vectors)."""
    phi = GOLDEN
    raw = []
    for s1 in (+1, -1):
        for s2 in (+1, -1):
            raw += [(0, s1*1, s2*phi), (s1*1, s2*phi, 0), (s1*phi, 0, s2*1)]
    v = np.unique(np.round(raw, 12), axis=0)
    v = v / np.linalg.norm(v, axis=1, keepdims=True)
    assert v.shape == (12, 3), v.shape
    return v

def icosahedral_selection(lmax=12):
    """Direct-a_lm power spectrum of 12 icosahedrally-placed delta sources. Returns P(l) normalised
    to its max. Theorem 1: P(l)~0 for l=1..5, O(1) at l=6."""
    v = icosahedron_vertices()
    theta = np.arccos(np.clip(v[:, 2], -1, 1))
    phi = np.arctan2(v[:, 1], v[:, 0])
    P = np.zeros(lmax + 1)
    for l in range(lmax + 1):
        s = 0.0
        for m in range(-l, l + 1):
            alm = np.sum(_Ylm(l, m, theta, phi))     # sum of delta sources -> a_lm = sum_j Ylm*(n_j)
            s += abs(alm) ** 2
        P[l] = s
    return P / P.max()


# ================================================================ Theorem 2
def ck_symmetric_alm(k, nside=128, lmax=8, seed=0):
    """a_lm of a C_k-symmetric field about z: sum of an off-axis blob over k azimuthal copies."""
    rng = np.random.default_rng(seed)
    npix = hp.nside2npix(nside)
    th, ph = hp.pix2ang(nside, np.arange(npix))
    # a few random blobs, each replicated at k azimuths -> exact C_k symmetry (up to pixelization)
    m = np.zeros(npix)
    for _ in range(4):
        th0 = rng.uniform(0.3, np.pi - 0.3); ph0 = rng.uniform(0, 2*np.pi); w = rng.uniform(0.15, 0.35)
        for j in range(k):
            pj = ph0 + 2*np.pi*j/k
            ang = np.arccos(np.clip(np.cos(th)*np.cos(th0)
                                    + np.sin(th)*np.sin(th0)*np.cos(ph - pj), -1, 1))
            m += np.exp(-0.5*(ang/w)**2)
    alm = hp.map2alm(m, lmax=lmax, iter=3)
    return alm

def m_power_by_residue(alm, k, lmax=8):
    """Return dict: residue r (m mod k) -> summed |a_lm|^2 over all l,m with that residue."""
    out = {r: 0.0 for r in range(k)}
    for l in range(lmax + 1):
        for mm in range(l + 1):
            out[mm % k] += abs(alm[hp.Alm.getidx(lmax, l, mm)]) ** 2
    tot = sum(out.values()) + 1e-300
    return {r: out[r] / tot for r in out}

def allowed_m_per_l(alm, lmax_show=4, tol=1e-3):
    """For each l, list m with fractional power above tol (in the z frame)."""
    res = {}
    for l in range(2, lmax_show + 1):
        pl = np.array([abs(alm[hp.Alm.getidx(hp.Alm.getlmax(alm.size), l, mm)])**2
                       for mm in range(l + 1)])
        pl = pl / (pl.sum() + 1e-300)
        res[l] = [int(mm) for mm in range(l + 1) if pl[mm] > tol]
    return res


# ================================================================ estimators
def _pure_Ylm_map(l, m, nside=64):
    """Real map of a single (l,m) harmonic (real part convention via healpy alm)."""
    lmax = max(l, 6)
    a = np.zeros(hp.Alm.getsize(lmax), dtype=complex)
    a[hp.Alm.getidx(lmax, l, m)] = 1.0
    return hp.alm2map(a, nside, verbose=False) if "verbose" in hp.alm2map.__code__.co_varnames \
        else hp.alm2map(a, nside)

def _power_tensor(mp, nside=64):
    vec = np.array(hp.pix2vec(nside, np.arange(hp.nside2npix(nside))))
    w = mp * mp
    return np.einsum('p,ip,jp->ij', w, vec, vec)

def axis_sectoral(mp, nside=64):
    """Smallest-eigval eigenvector (paper's sectoral / axis-of-evil estimator)."""
    ev, evec = np.linalg.eigh(_power_tensor(mp, nside)); n = evec[:, 0]; return n/np.linalg.norm(n)

def axis_zonal(mp, nside=64):
    """Largest-eigval eigenvector (paper's zonal estimator, section 5)."""
    ev, evec = np.linalg.eigh(_power_tensor(mp, nside)); n = evec[:, 2]; return n/np.linalg.norm(n)

def ang_to_z(n):
    return float(np.degrees(np.arccos(np.clip(abs(np.dot(n, [0, 0, 1.0])), 0, 1))))


# ================================================================ report
def _fmt(x): return ("%.2e" % x) if x < 1e-3 else ("%.4f" % x)

def main():
    print("="*74); print("THEOREM 1 -- icosahedral selection (direct a_lm of 12 face-centre sources)")
    P = icosahedral_selection(lmax=12)
    for l in range(1, 13):
        tag = "  <-- lowest invariant" if l == 6 else ""
        print("   P(l=%2d)/Pmax = %s%s" % (l, _fmt(P[l]), tag))
    ok1 = (P[1:6].max() < 1e-12) and (P[6] > 0.1)
    print("   VERDICT:", "PASS (l=1..5 vanish <1e-12; l=6 O(1))" if ok1 else "FAIL")

    print("="*74); print("THEOREM 2 -- k-fold azimuthal selection: surviving m are multiples of k")
    ok2 = True
    for k in (5, 3, 2):
        alm = ck_symmetric_alm(k, lmax=8, seed=1)
        frac = m_power_by_residue(alm, k, lmax=8)
        nonzero_res = {r: frac[r] for r in frac if r != 0}
        leak = max(nonzero_res.values()) if nonzero_res else 0.0
        print("   C%d: power in m==0(mod %d) = %.4f ; max leak into m!=0(mod %d) = %s"
              % (k, k, frac[0], k, _fmt(leak)))
        ok2 &= (frac[0] > 0.98) and (leak < 1e-2)
    print("   l-specific corollaries (|m|<=l with m==0 mod k):")
    print("     5-fold: l=2,3,4 -> m={0} (zonal);  first non-zonal at l=5,m=5")
    print("     3-fold: quadrupole -> m={0}; octopole -> m in {0,3} (m=+-3 allowed)")
    print("     2-fold: quadrupole -> m in {0,2} (m=+-2 allowed); octopole -> m in {0,2}, NO m=+-3")
    print("   VERDICT:", "PASS" if ok2 else "FAIL")

    print("="*74); print("APPENDIX A -- estimator self-tests")
    rows = []
    for (l, m, kind) in [(2, 2, "sectoral Y22"), (3, 3, "sectoral Y33"),
                         (2, 0, "zonal Y20"), (3, 0, "zonal Y30")]:
        mp = _pure_Ylm_map(l, m)
        es, ez = ang_to_z(axis_sectoral(mp)), ang_to_z(axis_zonal(mp))
        rows.append((kind, es, ez))
        print("   %-12s : sectoral-est err=%5.1f deg | zonal-est err=%5.1f deg" % (kind, es, ez))
    # expectation: sectoral recovers Y22/Y33 (~0) & misses Y20/Y30 (~90); zonal the reverse
    ok3 = (rows[0][1] < 2 and rows[1][1] < 2 and rows[2][1] > 88 and rows[3][1] > 88 and
           rows[0][2] > 88 and rows[1][2] > 88 and rows[2][2] < 2 and rows[3][2] < 2)
    print("   VERDICT:", "PASS (sectoral<->zonal are perpendicular, per section 4)" if ok3 else "FAIL")

    print("="*74); print("APPENDIX A -- five-fold H5 axis contrast (m=5 injection)")
    import eim_pipeline as E
    lmax = 8
    a = np.zeros(hp.Alm.getsize(lmax), dtype=complex)
    a[hp.Alm.getidx(lmax, 5, 5)] = 1.0; a[hp.Alm.getidx(lmax, 6, 5)] = 0.7
    mp = hp.alm2map(a, E.NSIDE)
    Hz = E.fivefold_H(mp, [0, 0, 1.0]); Hx = E.fivefold_H(mp, [1.0, 0, 0])
    ratio = Hz / (Hx + 1e-30)
    print("   H5(about z, true axis) = %.4f ; H5(about x) = %.4f ; z-to-x ratio = %.2f  (paper ~5.5)"
          % (Hz, Hx, ratio))
    ok4 = ratio > 3.0
    print("   VERDICT:", "PASS (strongly axis-selective)" if ok4 else "FAIL")

    print("="*74)
    allok = ok1 and ok2 and ok3 and ok4
    print("OVERALL:", "ALL SELECTION-RULE / ESTIMATOR CHECKS PASS" if allok else "SOME CHECKS FAILED")
    return allok

if __name__ == "__main__":
    main()
