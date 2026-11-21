"""Confirm01_JAG s80 (cloud): draw the N set per CONTRACT_CN.md section 2.2.
Usage: python s80_draw.py DataHoldings.txt tiles_145.txt OUTDIR
Region box and EPSG:3035 transform as in Probe_RelCons01 s51 (25<=lat<=75, -35<=lon<=50).
Strata are assigned in the order Al, S, C, Nt; within each stratum tiles are sorted by name and permuted with
numpy default_rng([20261117, k]) where k = 0,1,2,3 for Al,S,C,Nt. Quotas: Al 12, S 20, C 20, Nt 20."""
import sys, json
import numpy as np, pandas as pd, pyproj
DH, T145, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
old = set(x.strip() for x in open(T145) if x.strip()); assert len(old) == 145
rows = []
with open(DH, errors="ignore") as f:
    f.readline()
    for line in f:
        p = line.split()
        if len(p) < 9: continue
        try: lat, lon, h = float(p[1]), float(p[2]), float(p[3])
        except ValueError: continue
        if lon > 180: lon -= 360
        if 25 <= lat <= 75 and -35 <= lon <= 50:
            rows.append((p[0], lat, lon, h, p[7], p[8]))
s = pd.DataFrame(rows, columns=["sta", "lat", "lon", "h", "dtbeg", "dtend"]).drop_duplicates("sta")
s = s[(s.dtbeg <= "2020-01-01") & (s.dtend >= "2023-12-31")].copy()
x, y = pyproj.Transformer.from_crs(4326, 3035, always_xy=True).transform(s.lon.values, s.lat.values)
s["e"] = np.floor(x / 1e5).astype(int); s["n"] = np.floor(y / 1e5).astype(int)
s = s[(s.e >= 0) & (s.n >= 0)]
s["tile"] = [f"E{e}N{n}" for e, n in zip(s.e, s.n)]
s = s[~s.tile.isin(old)]
g = s.groupby("tile").agg(n=("n", "first"), n_sta=("sta", "size"), hmax=("h", "max"), stas=("sta", lambda v: " ".join(sorted(v)))).reset_index()
def stratum(r):
    if r.hmax >= 1000 and r.n <= 28: return "Al"
    if r.n <= 26: return "S"
    if r.n <= 35: return "C"
    return "Nt"
g["stratum"] = g.apply(stratum, axis=1)
Q = {"Al": 12, "S": 20, "C": 20, "Nt": 20}
out = []
for k, st in enumerate(["Al", "S", "C", "Nt"]):
    t = g[g.stratum == st].sort_values("tile").reset_index(drop=True)
    perm = np.random.default_rng([20261117, k]).permutation(len(t))
    t = t.iloc[perm].reset_index(drop=True); t["rank"] = np.arange(1, len(t) + 1)
    t["selected"] = t["rank"] <= Q[st]
    out.append(t)
o = pd.concat(out)[["stratum", "rank", "tile", "selected", "n_sta", "hmax", "stas"]]
o.to_csv(f"{OUT}/N_draw_order.csv", index=False)
json.dump({st: int((o.stratum == st).sum()) for st in Q}, open(f"{OUT}/N_draw_counts.json", "w"))
print(o.groupby("stratum").size().to_dict())
