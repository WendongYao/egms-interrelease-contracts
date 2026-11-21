"""GHSL exposure weights for RSC: verbatim slice of Confirm01 s89 (TR, ghsl_path, exposure)."""
import os, glob, zipfile
import numpy as np, pyproj, rasterio, rasterio.windows
GH = os.environ.get("GHSL_DIR", "")
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
