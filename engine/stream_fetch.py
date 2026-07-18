# ===CELL:md:intro===
"""
# EIM–CMB — streamlined data fetch (download → extract T → downgrade → compress)

**What this does.** Fetches the official Planck 2018 (PR3) component-separated CMB maps
(Commander, NILC, SEVEM, SMICA) + the common temperature mask, keeps **only the temperature
field**, harmonic-downgrades each to **Nside=64** (the engine's resolution), converts to µK,
removes monopole+dipole, and writes the small analysis inputs the frozen-axis engine consumes:

    smica_nside64_uK.fits   nilc_nside64_uK.fits
    sevem_nside64_uK.fits   commander_nside64_uK.fits
    commonmask_nside64.fits   gal070mask_nside64.fits   (TWO real masks; the decision rule needs >=2)
    planck_lowl_alm.npz      (ℓ≤32 alm per map — the maximally compressed archival form)
    planck_provenance.txt

**Why a streamlined version.** Each source map is ~1 GB (IQU, Nside=2048). The archive only
serves full-resolution IQU maps — there is no low-Nside SMICA/NILC product, and FITS BINTABLE
stores I/Q/U interleaved row-wise so you cannot range-request just the temperature column. The
download is therefore irreducible. What we *can* fix is the disk footprint: process **one map at
a time**, read **only field 0 (T)**, and **delete each ~1 GB raw file immediately** after
downgrading. Peak disk stays ~1 raw map (~1 GB) instead of the ~4 GB you'd get holding all four.

**Engine contract.** `run_real.py` reads only the five `*_nside64_uK.fits`/`commonmask_nside64.fits`
files and internally band-limits to ℓ≤32 at Nside=64. The `.npz` is an archival/quick-look
artifact (full-sky bandlimited alm); the engine re-derives masked alm from the FITS, so the npz is
not a substitute for the maps.

**Where to run.** Needs Planck-archive egress. Run in a session/machine with full network access.
The preflight cell checks both egress and free disk before downloading anything.
"""

# ===CELL:code:config===
import os, sys, time, shutil, tempfile, urllib.request, urllib.error
import numpy as np
import healpy as hp

# --- resolution / units (match the hardened engine: NSIDE=64, LMAX=32) ---
NSIDE_OUT = 64
LMAX_OUT  = 3 * NSIDE_OUT      # 192 — clean harmonic downgrade; engine only USES ℓ≤32
LMAX_USE  = 32                 # bandlimit stored in the .npz
K_TO_UK   = 1.0e6              # COM_CMB maps are in K_CMB

# --- where to write outputs (defaults to current dir) ---
WORKDIR = os.path.abspath(os.environ.get("EIM_WORKDIR", "."))
os.makedirs(WORKDIR, exist_ok=True)

# --- disk safety: require this much free space before each ~1 GB download ---
MIN_FREE_GB = 3.0

# --- archive locations ---
IRSA = "https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/cmb/"
PLA  = "https://pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID="

# --- the four PR3 component-separated full-mission maps (confirmed R3.00 filenames) ---
MAPS = {
    "smica":     "COM_CMB_IQU-smica_2048_R3.00_full.fits",
    "nilc":      "COM_CMB_IQU-nilc_2048_R3.00_full.fits",
    "sevem":     "COM_CMB_IQU-sevem_2048_R3.00_full.fits",
    "commander": "COM_CMB_IQU-commander_2048_R3.00_full.fits",
}
MASK = "COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits"   # single field TMASK, f_sky~0.78

# --- SECOND real mask (decision rule needs >=2). Planck Galactic-plane mask, apodized 2deg,
#     8 columns GAL020..GAL099; we take GAL070 (field index 3, f_sky~0.70) as a genuinely
#     different sky cut from the common mask. Change MASK2_FIELD to 4 for GAL080 (f_sky~0.80).
#     It lives under the PR2 ancillary tree; we try several candidate URLs + PLA and SKIP
#     gracefully (non-fatal) if unreachable, so a path change never kills the run. ---
MASK2 = "HFI_Mask_GalPlane-apo2_2048_R2.00.fits"
MASK2_FIELD = 3                                          # 0:GAL020 1:GAL040 2:GAL060 3:GAL070 4:GAL080 ...
MASK2_LABEL = "gal070"
_A2 = "https://irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/"
MASK2_URLS = [_A2 + "masks/" + MASK2, _A2 + MASK2, PLA + MASK2]

print("WORKDIR         :", WORKDIR)
print("maps to fetch   :", ", ".join(MAPS))
print("target Nside    :", NSIDE_OUT, "| downgrade lmax:", LMAX_OUT, "| npz bandlimit ℓ≤", LMAX_USE)

# ===CELL:code:helpers===
def free_gb(path=WORKDIR):
    """Free space (GB) on the filesystem holding `path`."""
    return shutil.disk_usage(path).free / 1e9

def _human(nbytes):
    for u in ("B", "KB", "MB", "GB"):
        if nbytes < 1024 or u == "GB":
            return f"{nbytes:.1f} {u}"
        nbytes /= 1024

def stream_download(fname, min_bytes=10_000_000, chunk=1 << 20, urls=None):
    """
    Stream `fname` to WORKDIR with a progress line, trying each URL in `urls` in order
    (default: IRSA then PLA on the component-maps path). Downloads to a .part temp then
    atomically renames, so a killed download never leaves a truncated 'final' file. Skips
    if the final file already exists. Returns the local path.
    """
    dst = os.path.join(WORKDIR, fname)
    if os.path.exists(dst) and os.path.getsize(dst) > min_bytes:
        print(f"  have {fname} ({_human(os.path.getsize(dst))}) — skip")
        return dst
    if free_gb() < MIN_FREE_GB:
        raise RuntimeError(
            f"Only {free_gb():.1f} GB free (< {MIN_FREE_GB} GB). Free space or lower MIN_FREE_GB.")
    if urls is None:
        urls = [IRSA + fname, PLA + fname]
    last_err = None
    for base in urls:
        part = dst + ".part"
        try:
            print(f"  downloading {fname} from {base.split('//')[1].split('/')[0]} ...")
            t0 = time.time()
            with urllib.request.urlopen(base, timeout=60) as r:
                total = int(r.headers.get("Content-Length", 0))
                got = 0
                with open(part, "wb") as fh:
                    while True:
                        buf = r.read(chunk)
                        if not buf:
                            break
                        fh.write(buf); got += len(buf)
                        if total:
                            pct = 100 * got / total
                            sys.stdout.write(
                                f"\r    {_human(got)} / {_human(total)} ({pct:4.1f}%)")
                        else:
                            sys.stdout.write(f"\r    {_human(got)}")
                        sys.stdout.flush()
            sys.stdout.write("\n")
            if os.path.getsize(part) <= min_bytes:
                raise RuntimeError(f"downloaded file too small ({os.path.getsize(part)} B)")
            os.replace(part, dst)
            print(f"    done in {time.time()-t0:.0f}s -> {_human(os.path.getsize(dst))}")
            return dst
        except Exception as e:               # noqa: BLE001 — try next mirror
            last_err = e
            print(f"    failed: {e}")
            if os.path.exists(part):
                os.remove(part)
    raise RuntimeError(f"could not download {fname} from any candidate URL: {last_err}")

def harmonic_downgrade_T(fname, field=0):
    """
    Read ONLY the temperature field, remove monopole+dipole on the good sky, harmonic-downgrade
    to Nside=64, convert to µK. Returns (map_uK[Nside=64], alm_lmaxUSE_uK_compact).

    Numerics are identical to the original fetch_planck.py downgrade (do not alter — the engine's
    matched-null validity assumes this exact anisotropy-field construction).
    """
    m = hp.read_map(fname, field=field, dtype=np.float64)          # K_CMB, field 0 = I only
    m = np.where(np.isfinite(m) & (np.abs(m) < 1e10), m, hp.UNSEEN)
    good = m != hp.UNSEEN
    mono, dip = hp.fit_dipole(hp.ma(np.where(good, m, hp.UNSEEN)))
    x, y, z = hp.pix2vec(hp.get_nside(m), np.arange(m.size))
    mm = np.where(good, m - mono - (dip[0]*x + dip[1]*y + dip[2]*z), 0.0)
    alm = hp.map2alm(mm, lmax=LMAX_OUT, iter=3)
    out = hp.alm2map(alm, NSIDE_OUT) * K_TO_UK                     # µK map at Nside=64
    # compact ℓ≤LMAX_USE alm of the µK map (true minimal representation, ~9 KB)
    alm32 = hp.map2alm(out, lmax=LMAX_USE, iter=3)
    return out, alm32

def lcdm_sanity(mp):
    """D_ℓ (µK²) at ℓ=2,5,10,20 from the Nside=64 µK map — expect an O(200–1200) plateau."""
    cl = hp.anafast(mp, lmax=64)
    ell = np.arange(cl.size)
    Dl = ell * (ell + 1) * cl / (2 * np.pi)
    return {int(i): round(float(Dl[i]), 0) for i in (2, 5, 10, 20)}

# ===CELL:code:preflight===
def preflight():
    ok = True
    # --- network egress to the archive ---
    try:
        req = urllib.request.Request(IRSA, method="HEAD")
        code = urllib.request.urlopen(req, timeout=20).status
        print(f"egress  : IRSA reachable (HTTP {code})")
    except Exception as e:                    # noqa: BLE001
        print(f"egress  : IRSA NOT reachable ({e})")
        print("          -> run this where the Planck archive is reachable (full network egress).")
        ok = False
    # --- disk ---
    fg = free_gb()
    est = 1.1 * len(MAPS) + 1.0               # rough GB touched over the run (one raw at a time)
    print(f"disk    : {fg:.1f} GB free (need ≥ {MIN_FREE_GB} GB headroom per download;"
          f" ~{est:.0f} GB will pass through, but peak resident ≈ 1 raw map)")
    if fg < MIN_FREE_GB:
        print("          -> insufficient free space for even one raw map."); ok = False
    print("PREFLIGHT:", "PASS" if ok else "FAIL — resolve the above before running the fetch cell.")
    return ok

# ===CELL:code:main===
def run_fetch():
    lowl_alm = {}
    for k, fname in MAPS.items():
        print(f"\n=== {k.upper()} ===")
        raw = stream_download(fname)                       # ~1 GB
        try:
            mp, alm32 = harmonic_downgrade_T(raw)
            hp.write_map(os.path.join(WORKDIR, f"{k}_nside64_uK.fits"), mp, overwrite=True)
            lowl_alm[k] = alm32
            print(f"   Nside=64 µK written | D_ℓ sanity {lcdm_sanity(mp)}  (LCDM ~200–1200)")
        finally:
            if os.path.exists(raw):                        # delete the ~1 GB raw NOW
                os.remove(raw)
                print(f"   removed raw {os.path.basename(raw)}  | free now {free_gb():.1f} GB")

    print("\n=== MASK ===")
    rawm = stream_download(MASK)
    try:
        mask = hp.read_map(rawm, dtype=np.float64)
        mask64 = (hp.ud_grade(mask, NSIDE_OUT) > 0.5).astype(np.float64)
        hp.write_map(os.path.join(WORKDIR, "commonmask_nside64.fits"), mask64, overwrite=True)
        print(f"   commonmask_nside64.fits written | f_sky = {mask64.mean():.3f}")
    finally:
        if os.path.exists(rawm):
            os.remove(rawm)
            print(f"   removed raw {os.path.basename(rawm)} | free now {free_gb():.1f} GB")

    print("\n=== SECOND MASK (%s, %s) ===" % (MASK2_LABEL, MASK2))
    mask2_ok = False
    try:
        rawm2 = stream_download(MASK2, urls=MASK2_URLS)     # ~1.6 GB, Galactic-plane, 8 columns
        try:
            m2 = hp.read_map(rawm2, field=MASK2_FIELD, dtype=np.float64)   # GAL070 column only
            m2_64 = (hp.ud_grade(m2, NSIDE_OUT) > 0.5).astype(np.float64)
            hp.write_map(os.path.join(WORKDIR, f"{MASK2_LABEL}mask_nside64.fits"), m2_64, overwrite=True)
            print(f"   {MASK2_LABEL}mask_nside64.fits written | f_sky = {m2_64.mean():.3f}")
            mask2_ok = True
        finally:
            if os.path.exists(rawm2):
                os.remove(rawm2)
                print(f"   removed raw {os.path.basename(rawm2)} | free now {free_gb():.1f} GB")
    except Exception as e:                                   # non-fatal: closure driver will warn
        print(f"   WARNING: second mask unreachable ({e}).")
        print(f"   -> run_closure.py will fall back to a synthetic second mask (still >=2, but the")
        print(f"      second is not a real Planck mask). Fix MASK2_URLS or supply the file manually.")

    # compressed archival alm + provenance
    np.savez(os.path.join(WORKDIR, "planck_lowl_alm.npz"), **lowl_alm)
    with open(os.path.join(WORKDIR, "planck_provenance.txt"), "w") as fh:
        fh.write("Planck 2018 PR3 component-separated CMB maps + temperature masks.\n")
        fh.write("Maps: " + ", ".join(MAPS.values()) + "\n")
        fh.write("Mask 1 (common): " + MASK + "\n")
        fh.write("Mask 2 (%s): %s field %d%s\n"
                 % (MASK2_LABEL, MASK2, MASK2_FIELD, "" if mask2_ok else "  [NOT FETCHED]"))
        fh.write("Temperature field only; monopole+dipole removed on the good sky.\n")
        fh.write(f"Harmonic-downgraded to Nside={NSIDE_OUT}, lmax={LMAX_OUT}; units µK_CMB.\n")
        fh.write(f".npz holds compact ℓ≤{LMAX_USE} alm (µK) per map. Engine uses ℓ≤32.\n")
    print("\nwrote planck_lowl_alm.npz + planck_provenance.txt")
    return lowl_alm

# ===CELL:code:verify===
def verify():
    """Reload the small outputs and check them — the deliverable's self-test."""
    import glob
    print("VERIFICATION")
    fits_out = [f"{k}_nside64_uK.fits" for k in MAPS] + ["commonmask_nside64.fits"]
    second = f"{MASK2_LABEL}mask_nside64.fits"
    if os.path.exists(os.path.join(WORKDIR, second)):
        fits_out.append(second)                       # include the real second mask when present
    allok = True
    for f in fits_out:
        p = os.path.join(WORKDIR, f)
        if not os.path.exists(p):
            print(f"  MISSING {f}"); allok = False; continue
        m = hp.read_map(p, dtype=np.float64, verbose=False) if "verbose" in hp.read_map.__code__.co_varnames \
            else hp.read_map(p, dtype=np.float64)
        ns = hp.get_nside(m); finite = np.isfinite(m).all()
        extra = ""
        if "mask" in f:
            vals = set(np.unique(m).tolist()); extra = f" f_sky={m.mean():.3f} binary={vals<= {0.0,1.0}}"
        else:
            extra = f" D_ℓ={lcdm_sanity(m)}"
        okrow = (ns == NSIDE_OUT) and finite
        allok &= okrow
        print(f"  {'OK ' if okrow else 'BAD'} {f:28s} Nside={ns} finite={finite}{extra}")
    # npz
    pnpz = os.path.join(WORKDIR, "planck_lowl_alm.npz")
    if os.path.exists(pnpz):
        d = np.load(pnpz)
        exp = hp.Alm.getsize(LMAX_USE)
        for k in MAPS:
            if k in d:
                good = d[k].size == exp
                allok &= good
                print(f"  {'OK ' if good else 'BAD'} npz[{k}] alm size={d[k].size} (expect {exp} for ℓ≤{LMAX_USE})")
    total = sum(os.path.getsize(os.path.join(WORKDIR, f))
                for f in fits_out + ["planck_lowl_alm.npz", "planck_provenance.txt"]
                if os.path.exists(os.path.join(WORKDIR, f)))
    print(f"  total output footprint: {_human(total)}  (vs ~{1.1*len(MAPS):.0f} GB of source downloads)")
    print("VERDICT:", "ALL CHECKS PASS — ready for run_real.py" if allok else "FAILURES ABOVE — do not proceed")
    return allok
