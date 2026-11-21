"""Release Stability Check (RSC) - core functions (Confirm02_Mechanism, CONTRACT section 4 Part III).
Pure functions on aligned arrays; I/O lives in rsc.py. No ranking of releases, no prediction.
  align(dfA, dfB)            coordinate-key alignment (rint(easting)*1e7 + rint(northing); first occurrence), as Confirm01 s82
  harmonize(...)             C0 / C1 (subtract per-point dg) / C1c (subtract a constant) / C2 (C1 + tile plane) / C2' (plane only)
  plane(E, N, r)             least-squares plane a + b*E + c*N on finite points (coordinates centred, 100 km units)
  hotspots_grid(E, N, mov)   moving points on a 100 m grid anchored at the data extent, 8-connected, clusters < 5 points dropped
  decision_metrics(...)      A1 / A2 / A3 for T5, T10, T0, T05 (definitions identical to Confirm01 s89)
"""
import numpy as np
from scipy import ndimage
TH = ["T5", "T10", "T0", "T05"]
def key(E, N):
    return np.rint(np.asarray(E, float)).astype(np.int64) * 10**7 + np.rint(np.asarray(N, float)).astype(np.int64)
def align(EA, NA, EB, NB):
    """Index arrays (ia, ib) of the points present in both releases, ordered by release A's first occurrences."""
    kA, kB = key(EA, NA), key(EB, NB)
    uA, fA = np.unique(kA, return_index=True); uB, fB = np.unique(kB, return_index=True)
    order = np.sort(fA); kk = kA[order]
    j = np.searchsorted(uB, kk); j = np.clip(j, 0, len(uB) - 1); ok = uB[j] == kk
    return order[ok], fB[j[ok]]
def plane(E, N, r):
    E = np.asarray(E, float); N = np.asarray(N, float); r = np.asarray(r, float); ok = np.isfinite(r)
    if ok.sum() < 3: return np.zeros_like(r)
    e0, n0 = E[ok].mean(), N[ok].mean(); X = np.c_[np.ones(len(E)), (E - e0) / 1e5, (N - n0) / 1e5]
    b = np.linalg.lstsq(X[ok], r[ok], rcond=None)[0]
    return X @ b
def harmonize(E, N, vA, vB, dg=None, mode="C0", const=0.0):
    """Return the harmonized release-B velocities for comparison mode C0, C1, C1c, C2 or C2p (C2')."""
    vA = np.asarray(vA, float); vB = np.asarray(vB, float)
    if mode == "C0": return vB.copy()
    if mode == "C1": return vB - np.asarray(dg, float)
    if mode == "C1c": return vB - float(const)
    if mode == "C2":
        b1 = vB - np.asarray(dg, float); return b1 - plane(E, N, b1 - vA)
    if mode == "C2p": return vB - plane(E, N, vB - vA)
    raise ValueError(mode)
def status(v, th):
    if th == "T5": return np.abs(v) >= 5
    if th == "T10": return np.abs(v) >= 10
    if th == "T0": return v < 0
    return np.where(v < -0.5, -1, np.where(v > 0.5, 1, 0))
def flagged(sa, sb, th):
    return (sa != 0) | (sb != 0) if th == "T05" else (sa | sb)
def hotspots_grid(E, N, mov):
    E = np.asarray(E); N = np.asarray(N)
    e0, n0 = (int(E.min()) // 100) * 100, (int(N.min()) // 100) * 100
    col = ((E - e0) // 100).astype(int); row = ((N - n0) // 100).astype(int)
    g = np.zeros((row.max() + 1, col.max() + 1), bool); g[row[mov], col[mov]] = True
    lab, k = ndimage.label(g, structure=np.ones((3, 3), int))
    if k == 0: return lab, [], row, col
    sizes = ndimage.sum(g, lab, index=np.arange(1, k + 1)); keep = [i + 1 for i, s in enumerate(sizes) if s >= 5]
    return lab, keep, row, col
def decision_metrics(E, N, vA, vB, built=None, pop=None, hotspot_table=False):
    """A1/A2/A3 for one area (tile or AOI). Points with non-finite vA or vB are dropped (as s89)."""
    vA = np.asarray(vA, float); vB = np.asarray(vB, float); ok = np.isfinite(vA) & np.isfinite(vB)
    vA, vB = vA[ok], vB[ok]; E, N = np.asarray(E)[ok], np.asarray(N)[ok]
    wb = np.zeros(len(vA)) if built is None else np.asarray(built, float)[ok]; wp = np.zeros(len(vA)) if pop is None else np.asarray(pop, float)[ok]
    res, tabs = {}, {}
    for th in TH:
        sa, sb = status(vA, th), status(vB, th); dif = sa != sb; fl = flagged(sa, sb, th)
        r = dict(n=int(ok.sum()), dis=float(dif.mean()), n_dis=int(dif.sum()), n_flag=int(fl.sum()),
                 built_dis=float(wb[dif].sum()), built_flag=float(wb[fl].sum()), pop_dis=float(wp[dif].sum()), pop_flag=float(wp[fl].sum()))
        if th in ("T5", "T10"):
            r.update(n_movA=int(sa.sum()), n_movB=int(sb.sum()), loss=float((sa & ~sb).sum() / sa.sum()) if sa.sum() else np.nan,
                     gain=float((sb & ~sa).sum() / sb.sum()) if sb.sum() else np.nan)
            la, ka, _, _ = hotspots_grid(E, N, sa); lb, kb, _, _ = hotspots_grid(E, N, sb)
            setb = set(kb); persA, ious, rows = [], [], []
            for i in ka:
                mA = la == i; ov = set(np.unique(lb[mA]).tolist()) & setb; persA.append(bool(ov))
                if ov: ious.append(max(float((mA & (lb == j)).sum()) / float((mA | (lb == j)).sum()) for j in ov))
                rows.append(dict(release="A", label=int(i), cells=int(mA.sum()), status="both" if ov else "only_A"))
            seta = set(ka); persB = []
            for j in kb:
                p = bool(set(np.unique(la[lb == j]).tolist()) & seta); persB.append(p)
                rows.append(dict(release="B", label=int(j), cells=int((lb == j).sum()), status="both" if p else "only_B"))
            r.update(n_hotA=len(ka), n_hotB=len(kb), persA=float(np.mean(persA)) if ka else np.nan, persB=float(np.mean(persB)) if kb else np.nan,
                     iou_med=float(np.median(ious)) if ious else np.nan)
            tabs[th] = rows
        res[th] = r
    return (res, tabs) if hotspot_table else res
