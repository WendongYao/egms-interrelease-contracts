"""Confirm01_JAG s86 (user's machine VM, read-only on raw zips) = Probe_RelCons01 s52 generalised to a manifest
(Amendment 1, pair P2 or any pair): per tile, sample 250 coordinate-matched points (rng seed [2, index of the tile in
the sorted manifest]) and all matched points within 300 m of an NGL station; output format identical to s52.
Usage: python s86_sample.py MANIFEST.csv STATIONS_3035.csv OUTDIR [time_budget_s] [worker k/n]
Paths in the manifest are absolute or relative to ~/mnt ($MNT). Prints progress only (no statistics)."""
import glob, json, os, re, sys, time, zipfile
T0 = time.time()
import numpy as np, pandas as pd
MAN, STA, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
BUDGET = float(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] else 1e12
WK = tuple(int(x) for x in sys.argv[5].split("/")) if len(sys.argv) > 5 else (0, 1)
os.makedirs(OUT, exist_ok=True)
RADIUS = 300.0; NS = 250
MNT = os.environ.get("MNT", os.path.expanduser("~/mnt"))
man = pd.read_csv(MAN)
ab = lambda p: p if os.path.isabs(p) else os.path.join(MNT, p)
by = {r.tile: {"A": ab(r.pathA), "B": ab(r.pathB)} for r in man.itertuples()}
tiles = sorted(by)
order = [tiles[i] for i in np.random.default_rng(20261110).permutation(len(tiles))]
sta = pd.read_csv(STA)

def open_csv(fz):
    z = zipfile.ZipFile(fz); n = [x for x in z.namelist() if x.endswith(".csv")][0]
    head = pd.read_csv(z.open(n), nrows=0).columns
    return z, n, [c for c in head if re.fullmatch(r"\d{8}", c)]
def coords(fz):
    z, n, dcols = open_csv(fz)
    c = pd.read_csv(z.open(n), usecols=["easting", "northing"])
    e = c["easting"].to_numpy(np.float64); nn = c["northing"].to_numpy(np.float64)
    return np.rint(e).astype(np.int64) * 10**7 + np.rint(nn).astype(np.int64), e, nn, dcols
def rows(fz, idx, dcols):
    z, n, _ = open_csv(fz); want = np.sort(np.asarray(idx)); out = []; off = 0
    for ch in pd.read_csv(z.open(n), usecols=dcols, dtype={c: np.float32 for c in dcols}, chunksize=50000):
        m = want[(want >= off) & (want < off + len(ch))]
        if len(m):
            out.append(ch.iloc[m - off][dcols].to_numpy(np.float32))
        off += len(ch)
    return np.vstack(out) if out else np.zeros((0, len(dcols)), np.float32)

for pos, t in enumerate(order):
    if pos % WK[1] != WK[0]:
        continue
    if os.path.exists(f"{OUT}/{t}.json"):
        continue
    if time.time() - T0 > BUDGET:
        print("time budget reached; rerun to resume", flush=True); sys.exit(0)
    ti = tiles.index(t); t1 = time.time()
    ka, ea, na, da = coords(by[t]["A"]); kb, eb, nb, db = coords(by[t]["B"])
    meta = dict(tile=t, pos=pos, sorted_index=ti, seen=False, nA=int(len(ka)), nB=int(len(kb)),
                dupA=int(len(ka) - len(np.unique(ka))), dupB=int(len(kb) - len(np.unique(kb))), nA_dates=len(da), nB_dates=len(db))
    # coordinate-matched points (first occurrence in each release)
    ua, ia = np.unique(ka, return_index=True); ub, ib = np.unique(kb, return_index=True)
    common, ca, cb = np.intersect1d(ua, ub, return_indices=True)
    rowA, rowB = ia[ca], ib[cb]  # row positions of the common points in A and B, aligned
    meta["n_matched"] = int(len(common))
    if len(common) < 50:
        meta["excluded"] = "fewer than 50 matched points"
        json.dump(meta, open(f"{OUT}/{t}.json", "w")); print(pos, t, "excluded", flush=True); continue
    rng = np.random.default_rng([2, ti])
    pick = np.sort(rng.choice(len(common), size=min(NS, len(common)), replace=False))
    # near-station matched points
    st = sta[sta.tile == t]; near_idx, near_sta, near_dist = [], [], []
    ce, cn = ea[rowA], na[rowA]
    for r in st.itertuples():
        d = np.hypot(ce - r.E, cn - r.N); j = np.where(d <= RADIUS)[0]
        near_idx += j.tolist(); near_sta += [r.sta] * len(j); near_dist += d[j].tolist()
    near_idx = np.array(near_idx, dtype=np.int64)
    allc = np.unique(np.concatenate([pick, near_idx])) if len(near_idx) else pick
    A = rows(by[t]["A"], rowA[allc], da); B = rows(by[t]["B"], rowB[allc], db)
    # rows() returns in sorted row order; map back to the order of allc
    oa = np.argsort(np.argsort(rowA[allc])); ob = np.argsort(np.argsort(rowB[allc]))
    A = A[oa]; B = B[ob]
    posmap = {int(c): k for k, c in enumerate(allc)}
    rsel = np.array([posmap[int(c)] for c in pick]); nsel = np.array([posmap[int(c)] for c in near_idx], dtype=np.int64)
    np.savez_compressed(f"{OUT}/{t}.npz", dA=np.array(da), dB=np.array(db), A_rand=A[rsel], B_rand=B[rsel],
                        E_rand=ce[pick], N_rand=cn[pick], A_sta=A[nsel] if len(nsel) else np.zeros((0, len(da)), np.float32),
                        B_sta=B[nsel] if len(nsel) else np.zeros((0, len(db)), np.float32), sta=np.array(near_sta), sta_dist=np.array(near_dist),
                        E_sta=ce[near_idx] if len(near_idx) else np.zeros(0), N_sta=cn[near_idx] if len(near_idx) else np.zeros(0))
    meta.update(n_rand=int(len(pick)), n_stations_in_tile=int(len(st)), n_near_points=int(len(near_idx)),
                n_stations_with_points=int(len(set(near_sta))), e_min=float(ea.min()), e_max=float(ea.max()), n_min=float(na.min()), n_max=float(na.max()),
                seconds=round(time.time() - t1, 1))
    json.dump(meta, open(f"{OUT}/{t}.json", "w"))
    print(pos, t, "ok", meta["n_matched"], meta["n_near_points"], meta["seconds"], "s", flush=True)
print("pass complete", flush=True)
