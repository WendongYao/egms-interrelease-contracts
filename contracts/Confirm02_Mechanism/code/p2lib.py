"""Confirm02_Mechanism p2lib: the Part II functions of Confirm01 s89 (status, flagged, hotspots, exposure, tile_metrics),
copied verbatim by slicing s89_part2.py (slice SHA-256 recorded in RUN_LOG); only the command-line parsing is replaced by
the GHSL_DIR environment variable."""
import sys, os, glob, json, re, zipfile
import numpy as np, pandas as pd, pyproj, rasterio, rasterio.windows
from scipy import ndimage
GH = os.environ.get("GHSL_DIR", "")
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
