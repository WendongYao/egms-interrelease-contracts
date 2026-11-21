"""Confirm01_JAG s89 (cloud): Part II (A) - stability of screening decisions across releases, on the data sets
O (old tiles x P1), O' (old tiles x P2) and N' (new tiles x P2); CONTRACT section 4 Part II and Amendment 1 section 2.3.
Inputs per set: full-tile outputs of s82 (all coordinate-matched points). Tile inclusion as in s54: >= 60 epochs in W
for each release, W >= 1095 days, >= 50 matched points.
Statuses (each velocity variant, each release): T5 |v|>=5; T10 |v|>=10; T0 v<0; T05 three classes (v<-0.5, |v|<=0.5, v>0.5).
Velocity variants: pub (published mean_velocity, each release's own window) and cw (common-window primary model).
Readouts:
  A1 point-level disagreement (statuses differ), plus for T5/T10 loss = P(not moving in B | moving in A) and gain = P(not moving in A | moving in B).
  A2 hotspots (T5, T10): moving points on the tile's 100 m grid, 8-connected, clusters with < 5 points dropped; persistence_A =
     share of A clusters sharing >= 1 cell with a B cluster, persistence_B the reverse; IoU median over A clusters (best overlapping B cluster).
  A3 exposure-weighted disagreement: GHSL R2023A 2020 built-up surface (m2) and population of the 100 m Mollweide cell holding each point;
     sum of weights of points whose status differs / sum of weights of points flagged in either release (T5/T10: moving; T0: subsiding;
     T05: not stable).
  A4 disagreement curve for |v| >= T, T = 0..15 step 0.5.
Aggregation (not fixed by the contract; fixed here before any result is read): A1, A2 primary = median over tiles of tile-level values
(each tile equal weight), pooled-over-points values reported as secondary; A3 primary = pooled ratio of sums over all included tiles.
Agreement between sets (pairs O-O', O'-N'; and O-N' reported): |diff| <= 0.05 or ratio in [2/3, 3/2]. Directional expectation:
cw disagreement < pub disagreement for T5 and T10 (A1 primary) in each of the three sets.
Usage: python s89_part2.py OUT_DIR GHSL_DIR SETNAME=FULLTILE_DIR [SETNAME=FULLTILE_DIR ...]"""
import sys, os, glob, json, re, zipfile
import numpy as np, pandas as pd, pyproj, rasterio, rasterio.windows
from scipy import ndimage
OUT, GH = sys.argv[1], sys.argv[2]; SETS = dict(a.split("=", 1) for a in sys.argv[3:]); os.makedirs(OUT, exist_ok=True)
TH = ["T5", "T10", "T0", "T05"]; VAR = ["pub", "cw"]
TR = pyproj.Transformer.from_crs(3035, "ESRI:54009", always_xy=True)
_tif = {}
def ghsl_path(kind, r, c):
    """Path of the extracted GHSL GeoTIFF (extracted once); None if the tile is not available."""
    key = (kind, r, c)
    if key not in _tif:
        pat = {"built": "GHS_BUILT_S_E2020", "pop": "GHS_POP_E2020"}[kind]
        fz = glob.glob(f"{GH}/{pat}*_R{r}_C{c}.zip")
        if not fz: _tif[key] = None
        else:
            z = zipfile.ZipFile(fz[0]); tif = [n for n in z.namelist() if n.endswith(".tif")][0]; p = f"{GH}/_x/{tif}"
            if not os.path.exists(p): z.extract(tif, f"{GH}/_x")
            _tif[key] = p
    return _tif[key]
def exposure(E, N):
    """Built-up surface and population of the GHSL 100 m cell holding each point (windowed reads; memory-bounded).
    Implementation change after the analysis-script freeze (whole-raster caching ran out of memory); values identical."""
    x, y = TR.transform(E.astype(float), N.astype(float))
    r = ((9000000 - y) // 1e6).astype(int) + 1; c = ((x + 18041000) // 1e6).astype(int) + 1
    out = {"built": np.zeros(len(E)), "pop": np.zeros(len(E))}
    for (rr, cc) in set(zip(r.tolist(), c.tolist())):
        m = (r == rr) & (c == cc)
        for kind in out:
            p = ghsl_path(kind, rr, cc)
            if p is None: continue
            with rasterio.open(p) as src:
                tf = src.transform; col = ((x[m] - tf.c) / tf.a).astype(int); row = ((y[m] - tf.f) / tf.e).astype(int)
                ok = (row >= 0) & (row < src.height) & (col >= 0) & (col < src.width); v = np.zeros(m.sum())
                if ok.any():
                    r0, r1, c0, c1 = row[ok].min(), row[ok].max() + 1, col[ok].min(), col[ok].max() + 1
                    a = src.read(1, window=rasterio.windows.Window(c0, r0, c1 - c0, r1 - r0)).astype(np.float64); nd = src.nodata
                    if nd is not None: a[a == nd] = 0.0
                    a[~np.isfinite(a)] = 0.0; a[a < 0] = 0.0
                    v[ok] = a[row[ok] - r0, col[ok] - c0]
                out[kind][m] = v
    return out
def status(v, th):
    if th == "T5": return np.abs(v) >= 5
    if th == "T10": return np.abs(v) >= 10
    if th == "T0": return v < 0
    return np.where(v < -0.5, -1, np.where(v > 0.5, 1, 0))
def flagged(sa, sb, th):
    return (sa != 0) | (sb != 0) if th == "T05" else (sa | sb)
def hotspots(E, N, mov, tile):
    e0, n0 = [int(x) * 100000 for x in re.findall(r"\d+", tile)]
    col = np.clip(((E - e0) // 100).astype(int), 0, 999); row = np.clip(((N - n0) // 100).astype(int), 0, 999)
    g = np.zeros((1000, 1000), bool); g[row[mov], col[mov]] = True
    lab, k = ndimage.label(g, structure=np.ones((3, 3), int))
    if k == 0: return lab, []
    sizes = ndimage.sum(g, lab, index=np.arange(1, k + 1)); keep = [i + 1 for i, s in enumerate(sizes) if s >= 5]
    return lab, keep
def tile_metrics(z, tile, expo):
    res = {}
    for var in VAR:
        vA, vB = z[f"vA_{var}"].astype(float), z[f"vB_{var}"].astype(float); ok = np.isfinite(vA) & np.isfinite(vB)
        vA, vB = vA[ok], vB[ok]; E, N = z["E"][ok], z["N"][ok]; wb, wp = expo["built"][ok], expo["pop"][ok]
        for th in TH:
            sa, sb = status(vA, th), status(vB, th); dif = sa != sb; fl = flagged(sa, sb, th)
            r = dict(n=int(ok.sum()), dis=float(dif.mean()), n_dis=int(dif.sum()), n_flag=int(fl.sum()),
                     built_dis=float(wb[dif].sum()), built_flag=float(wb[fl].sum()), pop_dis=float(wp[dif].sum()), pop_flag=float(wp[fl].sum()))
            if th in ("T5", "T10"):
                r.update(n_movA=int(sa.sum()), n_movB=int(sb.sum()), loss=float((sa & ~sb).sum() / sa.sum()) if sa.sum() else np.nan,
                         gain=float((sb & ~sa).sum() / sb.sum()) if sb.sum() else np.nan)
                la, ka = hotspots(E, N, sa, tile); lb, kb = hotspots(E, N, sb, tile)
                setb = set(kb); persA, ious = [], []
                for i in ka:
                    mA = la == i; ov = set(np.unique(lb[mA]).tolist()) & setb; persA.append(bool(ov))
                    if ov: ious.append(max(float((mA & (lb == j)).sum()) / float((mA | (lb == j)).sum()) for j in ov))
                seta = set(ka); persB = [bool(set(np.unique(la[lb == j]).tolist()) & seta) for j in kb]
                r.update(n_hotA=len(ka), n_hotB=len(kb), persA=float(np.mean(persA)) if ka else np.nan, persB=float(np.mean(persB)) if kb else np.nan,
                         iou_med=float(np.median(ious)) if ious else np.nan)
            res[(var, th)] = r
        curve = {}
        for T in np.arange(0, 15.01, 0.5):
            curve[float(T)] = float(((np.abs(vA) >= T) != (np.abs(vB) >= T)).mean())
        res[(var, "curve")] = curve
    return res
summ = {}; per_set_tiles = {}
for name, d in SETS.items():
    rows = []
    for jf in sorted(glob.glob(f"{d}/*.json")):
        m = json.load(open(jf)); t = m["tile"]
        inc = ("excluded" not in m) and m["nA_epochs"] >= 60 and m["nB_epochs"] >= 60 and m["W_days"] >= 1095 and m["n_matched"] >= 50
        if not inc: rows.append(dict(tile=t, included=False)); continue
        z = np.load(jf.replace(".json", ".npz")); expo = exposure(z["E"], z["N"])
        res = tile_metrics(z, t, expo)
        row = dict(tile=t, included=True, n=min(m["nA_epochs"], m["nB_epochs"]), Nidx=int(re.findall(r"N(\d+)", t)[0]))
        for (var, th), r in res.items():
            if th == "curve":
                for T, v in r.items(): row[f"{var}_curve_{T}"] = v
            else:
                for k, v in r.items(): row[f"{var}_{th}_{k}"] = v
        rows.append(row)
    T_ = pd.DataFrame(rows); T_.to_csv(f"{OUT}/by_tile_{name}.csv", index=False); per_set_tiles[name] = T_
    I = T_[T_.included == True].copy()
    I["dclass"] = np.where(I.n >= 200, "dense", np.where(I.n <= 130, "sparse", "mid")); I["band"] = np.where(I.Nidx <= 26, "S", np.where(I.Nidx <= 35, "C", "Nt"))
    S = dict(n_tiles=int(len(I)))
    for var in VAR:
        for th in TH:
            p = f"{var}_{th}_"
            s = dict(A1_tile_median=float(I[p + "dis"].median()), A1_pooled=float(I[p + "n_dis"].sum() / I[p + "n"].sum()),
                     A1_by_class={k: float(v[p + "dis"].median()) for k, v in I.groupby("dclass")},
                     A1_by_band={k: float(v[p + "dis"].median()) for k, v in I.groupby("band")},
                     A3_built=float(I[p + "built_dis"].sum() / I[p + "built_flag"].sum()) if I[p + "built_flag"].sum() > 0 else None,
                     A3_pop=float(I[p + "pop_dis"].sum() / I[p + "pop_flag"].sum()) if I[p + "pop_flag"].sum() > 0 else None)
            if th in ("T5", "T10"):
                s.update(loss_tile_median=float(I[p + "loss"].median()), gain_tile_median=float(I[p + "gain"].median()),
                         loss_pooled=float(((I[p + "loss"] * I[p + "n_movA"]).sum()) / I[p + "n_movA"].sum()) if I[p + "n_movA"].sum() else None,
                         A2_persA_tile_median=float(I[p + "persA"].median()), A2_persB_tile_median=float(I[p + "persB"].median()),
                         A2_iou_tile_median=float(I[p + "iou_med"].median()), n_hotA=int(I[p + "n_hotA"].sum()), n_hotB=int(I[p + "n_hotB"].sum()),
                         n_tiles_with_hotspots=int((I[p + "n_hotA"] > 0).sum()))
            S[f"{var}_{th}"] = s
        S[f"{var}_curve_tile_median"] = {T: float(I[f"{var}_curve_{T}"].median()) for T in np.arange(0, 15.01, 0.5).tolist()}
    S["direction_cw_lt_pub"] = {th: bool(S[f"cw_{th}"]["A1_tile_median"] < S[f"pub_{th}"]["A1_tile_median"]) for th in ("T5", "T10")}
    summ[name] = S
def agree(a, b):
    if a is None or b is None or (isinstance(a, float) and np.isnan(a)) or (isinstance(b, float) and np.isnan(b)): return None
    return bool(abs(a - b) <= 0.05 or (a > 0 and b > 0 and 2 / 3 <= a / b <= 1.5))
keys = []
for var in VAR:
    for th in TH:
        keys += [(f"{var}_{th}", "A1_tile_median"), (f"{var}_{th}", "A3_built"), (f"{var}_{th}", "A3_pop")]
        if th in ("T5", "T10"): keys += [(f"{var}_{th}", k) for k in ("loss_tile_median", "gain_tile_median", "A2_persA_tile_median", "A2_persB_tile_median", "A2_iou_tile_median")]
names = list(SETS); ag = []
for blk, k in keys:
    row = dict(item=f"{blk}.{k}", **{n: summ[n][blk].get(k) for n in names})
    for a, b in (("O", "O'"), ("O'", "N'"), ("O", "N'")):
        if a in summ and b in summ: row[f"agree_{a}_{b}"] = agree(summ[a][blk].get(k), summ[b][blk].get(k))
    row["conclusion"] = (bool(row.get("agree_O_O'") and row.get("agree_O'_N'") and row.get("agree_O_N'")) if all(n in summ for n in ("O", "O'", "N'")) else None)
    ag.append(row)
pd.DataFrame(ag).to_csv(f"{OUT}/agreement.csv", index=False)
summ["direction_all_sets"] = {th: bool(all(summ[n]["direction_cw_lt_pub"][th] for n in names)) for th in ("T5", "T10")}
json.dump(summ, open(f"{OUT}/summary.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("written", OUT)
