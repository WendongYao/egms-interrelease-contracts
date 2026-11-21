"""Confirm01_JAG s81: frozen B0/B1/B2 coefficients on the O set (CONTRACT section 3), using the same formulas as
Probe_EGMSUnc01 s72 (q95half; mvs fallback to the median positive mean_velocity_std) but fitted once on all
random-sample points of the 136 O tiles (no block folds). Usage: python s81_frozen_params.py BY_POINT_CSV OUT_JSON"""
import sys, json, hashlib
import numpy as np, pandas as pd
src, out = sys.argv[1], sys.argv[2]
P = pd.read_csv(src); P = P[P.kind == "rand"].copy(); P["absd"] = P.d.abs()
def q95half(x): return float(np.quantile(x, 0.95) / 2)
mvs_med = float(np.nanmedian(P.mean_velocity_std.where(P.mean_velocity_std > 0)))
def mvs(s): return s.where(np.isfinite(s) & (s > 0), mvs_med)
res = dict(source=src, source_sha256=hashlib.sha256(open(src, "rb").read()).hexdigest(), n_points=int(len(P)),
           n_tiles=int(P.tile.nunique()), mvs_fallback=mvs_med,
           pooled=dict(B0=q95half(P.absd / P.se), B1=q95half(P.absd / mvs(P.mean_velocity_std)), B2=q95half(P.absd)), by_class={})
for c in ("dense", "mid", "sparse"):
    t = P[P.dclass == c]
    res["by_class"][c] = dict(n_points=int(len(t)), n_tiles=int(t.tile.nunique()), B0=q95half(t.absd / t.se),
                              B1=q95half(t.absd / mvs(t.mean_velocity_std)), B2=q95half(t.absd))
json.dump(res, open(out, "w"), indent=1)
print(json.dumps(res["by_class"], indent=1)); print("pooled", res["pooled"])
