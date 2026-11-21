"""Confirm03_Proximity s97 (VM): extract gnss_velocity_e of both releases for every full-tile point of the east-west
(E) files of O' and N' (CONTRACT section 5, D-1: "另读 gnss_velocity_e 得 Δg_E"). Identical to Confirm02 s94 except
for the column name: points are aligned to the s82 full-tile outputs of the E files with the s82 coordinate key
(rint(easting) * 10**7 + rint(northing); first occurrence of each coordinate). Output per tile: npz with gA, gB (float32,
same order as the full-tile npz) and a json with counts. Prints only counts; never prints g values or differences.
Resumable: tiles with an existing output are skipped.
Usage: python s97_extract_gE.py MNT_ROOT MANIFEST FULLTILE_DIR OUT_DIR [MAX_SECONDS]"""
import sys, os, json, time, zipfile
import numpy as np, pandas as pd
COL = "gnss_velocity_e"
ROOT, MAN, FD, OUT = sys.argv[1:5]; TMAX = float(sys.argv[5]) if len(sys.argv) > 5 else 1e9; os.makedirs(OUT, exist_ok=True); t0 = time.time()
def gcol(fz):
    z = zipfile.ZipFile(fz); n = [x for x in z.namelist() if x.endswith(".csv")][0]
    m = pd.read_csv(z.open(n), usecols=["easting", "northing", COL], dtype=np.float64)
    key = np.rint(m.easting.values).astype(np.int64) * 10**7 + np.rint(m.northing.values).astype(np.int64)
    u, first = np.unique(key, return_index=True)
    return u, m[COL].values[first]
done = 0
for r in pd.read_csv(MAN).itertuples():
    t = r.tile
    if os.path.exists(f"{OUT}/{t}.json"): continue
    if not os.path.exists(f"{FD}/{t}.json"): continue  # full-tile output not ready yet
    if time.time() - t0 > TMAX: break
    js = json.load(open(f"{FD}/{t}.json"))
    if "excluded" in js:
        json.dump(dict(tile=t, excluded=js["excluded"]), open(f"{OUT}/{t}.json", "w")); continue
    z = np.load(f"{FD}/{t}.npz"); k = z["E"].astype(np.int64) * 10**7 + z["N"].astype(np.int64)
    out = {}; meta = dict(tile=t, n_fulltile=int(len(k)))
    for side, p in (("A", r.pathA), ("B", r.pathB)):
        u, g = gcol(f"{ROOT}/{p}"); i = np.searchsorted(u, k); i = np.clip(i, 0, len(u) - 1); ok = u[i] == k
        out[f"g{side}"] = np.where(ok, g[i], np.nan).astype(np.float32); meta[f"n_found_{side}"] = int(ok.sum())
        meta[f"n_nan_{side}"] = int((~np.isfinite(out[f"g{side}"])).sum())
    meta["aligned"] = bool(meta["n_found_A"] == len(k) and meta["n_found_B"] == len(k))
    np.savez_compressed(f"{OUT}/{t}.npz", **out); json.dump(meta, open(f"{OUT}/{t}.json", "w")); done += 1
    print(t, meta["n_fulltile"], meta["n_found_A"], meta["n_found_B"], meta["aligned"], flush=True)
print("done_this_call", done, "seconds", round(time.time() - t0, 1))
