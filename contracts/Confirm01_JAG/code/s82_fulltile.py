"""Confirm01_JAG s82 (user's machine VM, read-only on raw zips): full-tile extraction for Parts II and III.
For every coordinate-matched point of a tile (first occurrence per coordinate, as in Probe_RelCons01 s52):
  published mean_velocity of A and B; common-window primary-model velocity and formal s.e. for A and B
  (rc_common.window / tyears / fit_primary definitions; vectorised for gap-free rows, per-row fallback otherwise);
  mean_velocity_std (A, B), height (A); B3 quantities on the common dates (lag-1 autocorrelation of the
  detrended+deseasoned B-A difference series; correlation of the A and B primary-model residuals).
Tile metadata: window W, epochs in W per release, longest gap within W (days) per release, point counts.
Usage: python s82_fulltile.py MANIFEST.csv OUTDIR [time_budget_s] [worker k/n]
MANIFEST.csv columns: tile,set,pathA,pathB (absolute, or relative to ~/mnt or $MNT). Resumable. Prints progress only (no statistics)."""
import json, math, os, re, sys, time, zipfile
T0 = time.time()
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rc_common import window, tyears, fit_primary, d2ord
MAN, OUT = sys.argv[1], sys.argv[2]
BUDGET = float(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] else 1e12
WK = tuple(int(x) for x in sys.argv[4].split("/")) if len(sys.argv) > 4 else (0, 1)
os.makedirs(OUT, exist_ok=True)
META = ["easting", "northing", "height", "mean_velocity", "mean_velocity_std"]
CH = 40000
BLK = 50000  # rows per block for fitting (memory cap; results identical to whole-tile fitting)

def open_csv(fz):
    z = zipfile.ZipFile(fz); n = [x for x in z.namelist() if x.endswith(".csv")][0]
    head = pd.read_csv(z.open(n), nrows=0).columns
    return z, n, [c for c in head if re.fullmatch(r"\d{8}", c)]

def read_meta(fz):
    z, n, dcols = open_csv(fz)
    head = list(pd.read_csv(z.open(n), nrows=0).columns); use = [c for c in META if c in head]
    m = pd.read_csv(z.open(n), usecols=use, dtype={c: np.float64 for c in use})
    for c in META:
        if c not in m: m[c] = np.nan  # e.g. some releases carry no height column
    key = np.rint(m.easting.values).astype(np.int64) * 10**7 + np.rint(m.northing.values).astype(np.int64)
    return m, key, dcols

def read_rows(fz, rows, cols):
    """Rows (sorted positions) x cols as float64, streaming in chunks; returns array aligned to `rows` order."""
    z, n, _ = open_csv(fz); order = np.argsort(rows); srt = rows[order]
    out = np.empty((len(rows), len(cols)), np.float32); off = 0; k = 0
    for ch in pd.read_csv(z.open(n), usecols=cols, dtype={c: np.float32 for c in cols}, chunksize=CH):
        lo = np.searchsorted(srt, off); hi = np.searchsorted(srt, off + len(ch))
        if hi > lo:
            out[order[lo:hi]] = ch[cols].to_numpy(np.float32)[srt[lo:hi] - off]
        off += len(ch)
    return out

def design(t): return np.vstack([np.ones_like(t), t, np.sin(2 * math.pi * t), np.cos(2 * math.pi * t)]).T

def fit_all(t, Y):
    """Primary model per row: returns v, se, residual matrix (NaN where not computable)."""
    n = Y.shape[0]; v = np.full(n, np.nan); se = np.full(n, np.nan); R = np.full_like(Y, np.nan)
    X = design(t); fin = np.isfinite(Y); full = fin.all(1)
    if full.any() and len(t) >= 60:
        P = np.linalg.pinv(X); C = Y[full] @ P.T; Rf = Y[full] - C @ X.T
        s2 = (Rf ** 2).sum(1) / max(len(t) - 4, 1); xtx = np.linalg.pinv(X.T @ X)[1, 1]
        v[full] = C[:, 1]; se[full] = np.sqrt(np.maximum(s2 * xtx, 0)); R[full] = Rf
    for i in np.where(~full)[0]:
        r = fit_primary(t, Y[i])
        if r is not None: v[i], se[i] = r
    return v, se, R

def max_gap(dates):
    o = np.array([d2ord(x) for x in dates]); return int(np.diff(o).max()) if len(o) > 1 else -1

man = pd.read_csv(MAN)
MNT = os.environ.get("MNT", os.path.expanduser("~/mnt"))
for c in ("pathA", "pathB"): man[c] = [p if os.path.isabs(p) else os.path.join(MNT, p) for p in man[c]]
for pos, r in man.iterrows():
    if pos % WK[1] != WK[0]: continue
    t = r.tile; fo = f"{OUT}/{t}.npz"
    if os.path.exists(fo) and os.path.exists(f"{OUT}/{t}.json"): continue
    if time.time() - T0 > BUDGET: print("time budget reached; rerun to resume", flush=True); sys.exit(0)
    t1 = time.time()
    mA, kA, dA = read_meta(r.pathA); mB, kB, dB = read_meta(r.pathB)
    ua, ia_ = np.unique(kA, return_index=True); ub, ib_ = np.unique(kB, return_index=True)
    _, ca, cb = np.intersect1d(ua, ub, return_indices=True); rowA, rowB = ia_[ca], ib_[cb]
    lo, hi, ia, ib = window(dA, dB)
    wA = [dA[i] for i in ia]; wB = [dB[i] for i in ib]
    meta = dict(tile=t, set=r.set, W_start=lo, W_end=hi, W_days=d2ord(hi) - d2ord(lo), nA_epochs=len(wA), nB_epochs=len(wB),
                n=min(len(wA), len(wB)), gapA_days=max_gap(wA), gapB_days=max_gap(wB), nA_points=int(len(kA)), nB_points=int(len(kB)),
                dupA=int(len(kA) - len(ua)), dupB=int(len(kB) - len(ub)), n_matched=int(len(rowA)), fileA=os.path.basename(r.pathA), fileB=os.path.basename(r.pathB))
    if len(rowA) == 0 or len(wA) < 4 or len(wB) < 4:
        meta["excluded"] = "no matched points or empty window"; json.dump(meta, open(f"{OUT}/{t}.json", "w")); continue
    YA = read_rows(r.pathA, rowA, wA); YB = read_rows(r.pathB, rowB, wB)
    common = sorted(set(wA) & set(wB)); jA = [wA.index(x) for x in common]; jB = [wB.index(x) for x in common]
    nr = len(rowA); vAc = np.full(nr, np.nan); seA = np.full(nr, np.nan); vBc = np.full(nr, np.nan); seB = np.full(nr, np.nan)
    ac1 = np.full(nr, np.nan); rc = np.full(nr, np.nan); tA, tB = tyears(wA), tyears(wB); tc = tyears(common) if common else None
    def corr_rows(a, b):
        a = a - a.mean(1, keepdims=True); b = b - b.mean(1, keepdims=True)
        return (a * b).sum(1) / np.sqrt((a * a).sum(1) * (b * b).sum(1))
    for b0 in range(0, nr, BLK):
        sl = slice(b0, min(nr, b0 + BLK)); ya = YA[sl].astype(np.float64); yb = YB[sl].astype(np.float64)
        vAc[sl], seA[sl], _ = fit_all(tA, ya); vBc[sl], seB[sl], _ = fit_all(tB, yb)
        if len(common) >= 60:
            D = yb[:, jB] - ya[:, jA]
            _, _, Rd = fit_all(tc, D); ok = np.isfinite(Rd).all(1); tmp = np.full(Rd.shape[0], np.nan)
            tmp[ok] = corr_rows(Rd[ok, :-1], Rd[ok, 1:]); ac1[sl] = tmp; del D, Rd
            _, _, Ra = fit_all(tc, ya[:, jA]); _, _, Rb = fit_all(tc, yb[:, jB]); ok2 = np.isfinite(Ra).all(1) & np.isfinite(Rb).all(1)
            tmp = np.full(Ra.shape[0], np.nan); tmp[ok2] = corr_rows(Ra[ok2], Rb[ok2]); rc[sl] = tmp; del Ra, Rb
        del ya, yb
    nfA = np.isfinite(YA).sum(1).astype(np.int16); nfB = np.isfinite(YB).sum(1).astype(np.int16); del YA, YB
    meta.update(n_common_dates=len(common), seconds=round(time.time() - t1, 1))
    f32 = lambda x: np.asarray(x, np.float32)
    np.savez_compressed(fo, E=np.rint(mA.easting.values[rowA]).astype(np.int32), N=np.rint(mA.northing.values[rowA]).astype(np.int32),
                        height=f32(mA.height.values[rowA]), vA_pub=f32(mA.mean_velocity.values[rowA]), vB_pub=f32(mB.mean_velocity.values[rowB]),
                        mvsA=f32(mA.mean_velocity_std.values[rowA]), mvsB=f32(mB.mean_velocity_std.values[rowB]),
                        vA_cw=f32(vAc), vB_cw=f32(vBc), seA=f32(seA), seB=f32(seB),
                        nfinA=nfA, nfinB=nfB, ac1_diff=f32(ac1), rescorr=f32(rc))
    json.dump(meta, open(f"{OUT}/{t}.json", "w"))
    print(f"{t}: matched={len(rowA)} epochs=({len(wA)},{len(wB)}) {meta['seconds']}s", flush=True)
print("done", flush=True)
