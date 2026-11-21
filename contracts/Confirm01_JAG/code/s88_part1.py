"""Confirm01_JAG s88 (cloud): Part I readouts F1-F11 for one data set (N' primary, O' control), CONTRACT section 4
and Amendment 1 section 2.2. Per-point and per-tile definitions copied from Probe_RelCons01 s54 (velocities, flips,
m, frac>2sigma, station medians, H2-type ratio) and Probe_EGMSUnc01 s72 (formal s.e. of release A, mean_velocity_std,
interval score, tercile ratio, spatial ratio T3, GNSS ordering T4). Coefficients are the frozen ones (FROZEN_PARAMS.json,
fitted on O x P1); nothing is refitted. Written and hashed before any extraction output of N' or O' is read.
Usage: python s88_part1.py SAMPLE_DIR FULLTILE_DIR GNSS_CSV FROZEN_JSON OUT_DIR"""
import sys, os, glob, json, math
import numpy as np, pandas as pd
from scipy.stats import spearmanr, binomtest
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rc_common import window, tyears, fit_ols, density_class, d2ord
from s31_realcf import year_fracs, rmse_lin_harm, labels, compare
SD, FD, GC, FZ, R = sys.argv[1:6]; os.makedirs(R, exist_ok=True)
SEED = 20261117; NBOOT = 1000
FR = json.load(open(FZ))

def fit_full(t, y, min_n=60):
    fin = np.isfinite(y)
    if fin.sum() < min_n: return None
    tt = t[fin]; yy = y[fin].astype(np.float64)
    X = np.vstack([np.ones_like(tt), tt, np.sin(2 * math.pi * tt), np.cos(2 * math.pi * tt)]).T
    coef, *_ = np.linalg.lstsq(X, yy, rcond=None); r = yy - X @ coef
    s2 = float(r @ r) / max(len(yy) - 4, 1); se = math.sqrt(max((s2 * np.linalg.pinv(X.T @ X))[1, 1], 0.0))
    return float(coef[1]), se

tile_rows, pt_rows, st_rows = [], [], []
for jf in sorted(glob.glob(f"{SD}/*.json")):
    m = json.load(open(jf)); t = m["tile"]
    row = dict(tile=t, n_matched=m.get("n_matched"), nA_points=m.get("nA"), nB_points=m.get("nB"))
    if "excluded" in m:
        row.update(included=False, reason=m["excluded"]); tile_rows.append(row); continue
    z = np.load(jf.replace(".json", ".npz"), allow_pickle=True)
    lo, hi, ia, ib = window(z["dA"], z["dB"]); dA = [str(z["dA"][i]) for i in ia]; dB = [str(z["dB"][i]) for i in ib]
    n = min(len(ia), len(ib)); wdays = d2ord(hi) - d2ord(lo) if (ia and ib) else 0
    row.update(W_start=lo, W_end=hi, W_days=wdays, nA_epochs=len(ia), nB_epochs=len(ib), n=n, dclass=density_class(n),
               northing_mean=float(np.mean(z["N_rand"])), easting_mean=float(np.mean(z["E_rand"])))
    if len(ia) < 60 or len(ib) < 60 or wdays < 1095:
        row.update(included=False, reason="window: <60 epochs or <3 years"); tile_rows.append(row); continue
    # mean_velocity_std (release A) from the full-tile output, matched on rounded coordinates
    ft = np.load(f"{FD}/{t}.npz"); key = ft["E"].astype(np.int64) * 10**7 + ft["N"].astype(np.int64)
    mvs_map = dict(zip(key.tolist(), ft["mvsA"].astype(float).tolist()))
    def mvs_of(E, N): return np.array([mvs_map.get(int(round(e)) * 10**7 + int(round(nn)), np.nan) for e, nn in zip(E, N)])
    tA, tB = tyears(dA), tyears(dB); yfA = year_fracs(dA)
    A, B = z["A_rand"][:, ia].astype(np.float64), z["B_rand"][:, ib].astype(np.float64); mv = mvs_of(z["E_rand"], z["N_rand"])
    recs = []
    for p in range(A.shape[0]):
        fa, fb = fit_full(tA, A[p]), fit_full(tB, B[p])
        if fa is None or fb is None: continue
        oa, ob = fit_ols(tA, A[p]), fit_ols(tB, B[p])
        fin = np.isfinite(A[p]); rmse = rmse_lin_harm(yfA[fin], A[p][fin])
        LA, LB = labels(A[p], dA, rmse), labels(B[p], dB, rmse)
        rec = dict(tile=t, kind="rand", sta="", point=p, E=float(z["E_rand"][p]), N=float(z["N_rand"][p]), vA=fa[0], vB=fb[0], seA=fa[1], seB=fb[1],
                   vA_ols=oa, vB_ols=ob, mean_velocity_std=float(mv[p]), dclass=density_class(n))
        if LA is not None and LB is not None:
            c = compare(LA, LB); rec.update(trend_changed=bool(c["A_class"]["changed"]), trend_signflip=bool(c["A_class"]["sign_flip"]),
                                            step_changed=bool(c["C_class"]["changed"]), seas_changed=bool(c["B_class"]["changed"]))
        recs.append(rec)
    P = pd.DataFrame(recs)
    for c_ in ("trend_changed", "trend_signflip", "step_changed", "seas_changed"):
        P[c_] = P[c_].astype(float) if c_ in P else np.nan
    if len(P) < 50:
        row.update(included=False, reason="fewer than 50 usable points"); tile_rows.append(row); continue
    pt_rows += P.to_dict("records")
    dv = P.vB - P.vA; dvo = P.vB_ols - P.vA_ols
    row.update(included=True, reason="", n_points=len(P), m=float(dv.abs().median()), m_ols=float(dvo.abs().median()), bias=float(dv.median()),
               p90_absdv=float(dv.abs().quantile(0.9)), frac_gt_2sigma=float((dv.abs() > 2 * np.sqrt(P.seA ** 2 + P.seB ** 2)).mean()),
               trend_changed=float(P.trend_changed.mean()), step_changed=float(P.step_changed.mean()), seas_changed=float(P.seas_changed.mean()),
               n_lab=int(P.trend_changed.notna().sum()), n_trend_changed=int(P.trend_changed.sum()), n_step_changed=int(P.step_changed.sum()),
               n_seas_changed=int(P.seas_changed.sum()))
    tile_rows.append(row)
    if len(z["sta"]):
        As, Bs = z["A_sta"][:, ia].astype(np.float64), z["B_sta"][:, ib].astype(np.float64); names = np.array([str(s) for s in z["sta"]])
        mvs_s = mvs_of(z["E_sta"], z["N_sta"])
        for s in sorted(set(names)):
            va, vb = [], []
            for k in np.where(names == s)[0]:
                fa, fb = fit_full(tA, As[k]), fit_full(tB, Bs[k])
                if fa is not None and fb is not None:
                    va.append(fa[0]); vb.append(fb[0])
                    pt_rows.append(dict(tile=t, kind="sta", sta=s, E=float(z["E_sta"][k]), N=float(z["N_sta"][k]), vA=fa[0], vB=fb[0], seA=fa[1], seB=fb[1],
                                        mean_velocity_std=float(mvs_s[k]), dclass=density_class(n)))
            if len(va) >= 3:
                st_rows.append(dict(sta=s, tile=t, n=n, dclass=density_class(n), n_pts=len(va), vA=float(np.median(va)), vB=float(np.median(vb))))
T = pd.DataFrame(tile_rows); T.to_csv(f"{R}/by_tile.csv", index=False)
PT = pd.DataFrame(pt_rows); PT["d"] = PT.vB - PT.vA; PT["absd"] = PT.d.abs(); PT.to_csv(f"{R}/by_point.csv", index=False)
inc = T[T.included == True].copy()
S = dict(n_tiles=int(len(T)), n_included=int(len(inc)), excluded=T[T.included != True][["tile", "reason"]].to_dict("records"),
         by_class_tiles={k: int(v) for k, v in inc.dclass.value_counts().items()})
rng = lambda: np.random.default_rng(SEED)
# F1
nl = inc.n_lab.sum()
S["F1"] = dict(n_points=int(nl), trend=float(inc.n_trend_changed.sum() / nl), step=float(inc.n_step_changed.sum() / nl), seasonal=float(inc.n_seas_changed.sum() / nl))
S["F1"]["pass"] = bool(0.08 <= S["F1"]["trend"] <= 0.20 and 0.15 <= S["F1"]["step"] <= 0.35)
# F2
rho = spearmanr(inc.n, inc.m).statistic; g = rng(); bs = []
for _ in range(NBOOT):
    s = inc.iloc[g.integers(0, len(inc), len(inc))]; bs.append(spearmanr(s.n, s.m).statistic if s.n.nunique() > 1 else np.nan)
md, ms = inc[inc.dclass == "dense"].m, inc[inc.dclass == "sparse"].m; ratio_ok = len(md) >= 10 and len(ms) >= 10
S["F2"] = dict(rho=float(rho), rho_ci_lo=float(np.nanpercentile(bs, 2.5)), rho_ci_hi=float(np.nanpercentile(bs, 97.5)),
               med_by_class={k: float(v.median()) for k, v in inc.groupby("dclass").m}, ratio=float(ms.median() / md.median()) if (len(md) and len(ms)) else None,
               ratio_evaluated=bool(ratio_ok), rho_ols=float(spearmanr(inc.n, inc.m_ols).statistic))
S["F2"]["pass"] = bool(rho <= -0.5 and S["F2"]["rho_ci_hi"] < -0.2 and ((S["F2"]["ratio"] >= 3) if ratio_ok else True))
# F3
f3 = {k: float(v.median()) for k, v in inc.groupby("dclass").frac_gt_2sigma}
S["F3"] = dict(frac_gt_2sigma_by_class=f3, **{"pass": bool(all(v >= 0.5 for v in f3.values()))})
# F4 / F5 (random points of included tiles)
E_ = PT[(PT.kind == "rand") & PT.tile.isin(inc.tile)].copy()
mvs_med = FR["mvs_fallback"]
def mvs(s): return s.where(np.isfinite(s) & (s > 0), mvs_med)
k0 = E_.dclass.map({c: FR["by_class"][c]["B0"] for c in FR["by_class"]}); k1 = E_.dclass.map({c: FR["by_class"][c]["B1"] for c in FR["by_class"]})
E_["sig_B0"] = k0 * E_.seA; E_["sig_B1"] = k1 * mvs(E_.mean_velocity_std)
cov = dict(pooled=float((E_.absd <= 2 * E_.sig_B0).mean()), by_class={k: dict(n_tiles=int(g_.tile.nunique()), coverage=float((g_.absd <= 2 * g_.sig_B0).mean())) for k, g_ in E_.groupby("dclass")})
S["F4"] = dict(coverage=cov, **{"pass": bool(0.92 <= cov["pooled"] <= 0.97 and all(v["coverage"] >= 0.88 for v in cov["by_class"].values() if v["n_tiles"] >= 5))})
def iscore(sig, absd): return 4 * sig + 40 * np.maximum(0, absd - 2 * sig)
def tercile_ratio(g_, col):
    if g_[col].nunique() < 3: return np.nan
    lo_, hi_ = g_[col].quantile(1 / 3), g_[col].quantile(2 / 3); a, bb = g_[g_[col] <= lo_].absd.median(), g_[g_[col] >= hi_].absd.median()
    return bb / a if a > 0 else np.nan
is0, is1 = float(iscore(E_.sig_B0, E_.absd).mean()), float(iscore(E_.sig_B1, E_.absd).mean())
tr0 = float(np.nanmedian(E_.groupby("tile").apply(lambda g_: tercile_ratio(g_, "sig_B0")))); tr1 = float(np.nanmedian(E_.groupby("tile").apply(lambda g_: tercile_ratio(g_, "sig_B1"))))
ndist = int(E_.mean_velocity_std.dropna().round(1).nunique())
S["F5"] = dict(n_distinct_mvs=ndist, n_mvs_missing=int(E_.mean_velocity_std.isna().sum()), interval_score_B0=is0, interval_score_B1=is1, ratio_B1_over_B0=is1 / is0,
               tercile_ratio_B0=tr0, tercile_ratio_B1=tr1, coverage_B1=float((E_.absd <= 2 * E_.sig_B1).mean()))
S["F5"]["pass"] = bool(ndist <= 20 and is1 >= 1.10 * is0 and tr1 < tr0)
# F6 (spatial ratio, station clusters; T3 definition)
sp = []
for t, g_ in PT[(PT.kind == "sta") & PT.tile.isin(inc.tile)].groupby("tile"):
    cl = [(s, h) for s, h in g_.groupby("sta") if len(h) >= 3]
    if len(cl) < 2: continue
    cen = {s: (h.E.mean(), h.N.mean()) for s, h in cl}; short, long_ = [], []
    for i, (s1, h1) in enumerate(cl):
        d1 = h1.d.values; iu = np.triu_indices(len(d1), 1); short += np.abs(d1[:, None] - d1[None, :])[iu].tolist()
        for s2, h2 in cl[i + 1:]:
            if math.hypot(cen[s1][0] - cen[s2][0], cen[s1][1] - cen[s2][1]) > 20000:
                long_ += np.abs(d1[:, None] - h2.d.values[None, :]).ravel().tolist()
    if short and long_ and np.median(long_) > 0:
        sp.append(dict(tile=t, dclass=g_.dclass.iloc[0], n_clusters=len(cl), short=float(np.median(short)), long=float(np.median(long_)), r=float(np.median(short) / np.median(long_))))
SP = pd.DataFrame(sp); SP.to_csv(f"{R}/by_tile_spatial.csv", index=False)
if len(SP) >= 15:
    g = rng(); bs = [float(np.median(SP.r.values[g.integers(0, len(SP), len(SP))])) for _ in range(NBOOT)]
    S["F6"] = dict(n_tiles=int(len(SP)), r_median=float(SP.r.median()), r_ci_lo=float(np.percentile(bs, 2.5)), r_ci_hi=float(np.percentile(bs, 97.5)),
                   by_class={k: dict(n=int(len(v)), r_median=float(v.r.median())) for k, v in SP.groupby("dclass")}, frac_below_1=float((SP.r < 1).mean()))
    S["F6"]["pass"] = bool(S["F6"]["r_median"] <= 0.8 and S["F6"]["r_ci_hi"] < 1.0)
else:
    S["F6"] = dict(n_tiles=int(len(SP)), untestable=True, **{"pass": False})
# GNSS: F7-F10
ST = pd.DataFrame(st_rows)
if len(ST) and os.path.exists(GC):
    G = pd.read_csv(GC); G = G[G.included == True][["sta", "v_gnss", "v_gnss_se_white", "completeness", "n_steps"]]
    ST = ST.merge(G, on="sta", how="inner")
else:
    ST = ST.iloc[0:0]
nst = int(len(ST)); nsp = int((ST.dclass == "sparse").sum()) if nst else 0
S["n_stations"] = nst; S["n_stations_sparse"] = nsp
if nst >= 5:
    c0 = float(np.median((ST.vA + ST.vB) / 2 - ST.v_gnss)); ST["eA"] = np.abs(ST.vA - c0 - ST.v_gnss); ST["eB"] = np.abs(ST.vB - c0 - ST.v_gnss)
    ST["m_err"] = np.abs((ST.vA + ST.vB) / 2 - c0 - ST.v_gnss); ST["absdv"] = (ST.vB - ST.vA).abs()
    rho7 = float(spearmanr(ST.absdv, ST.m_err).statistic); g = rng(); bs = []
    for _ in range(NBOOT):
        s = ST.iloc[g.integers(0, nst, nst)]; bs.append(spearmanr(s.absdv, s.m_err).statistic)
    lo_, hi_ = ST.absdv.quantile(1 / 3), ST.absdv.quantile(2 / 3)
    S["F7"] = dict(rho=rho7, rho_p05=float(np.nanpercentile(bs, 5)), tercile_ratio=float(ST[ST.absdv >= hi_].m_err.median() / ST[ST.absdv <= lo_].m_err.median()),
                   by_class={k: float(spearmanr(v.absdv, v.m_err).statistic) for k, v in ST.groupby("dclass") if len(v) >= 10})
    S["F7"]["testable"] = bool(nst >= 60 and nsp >= 15); S["F7"]["pass"] = bool(S["F7"]["testable"] and rho7 >= 0.10 and S["F7"]["rho_p05"] > 0)
    def rstat(d):
        c = float(np.median((d.vA + d.vB) / 2 - d.v_gnss)); eA = np.abs(d.vA - c - d.v_gnss).values; eB = np.abs(d.vB - c - d.v_gnss).values
        den = float(np.median(np.concatenate([eA, eB]))); return float(np.median(np.abs(d.vB - d.vA))) / den if den > 0 else np.nan
    g = rng(); bs = [rstat(ST.iloc[g.integers(0, nst, nst)]) for _ in range(NBOOT)]
    S["F8"] = dict(R=rstat(ST), R_p05=float(np.nanpercentile(bs, 5)), R_p95=float(np.nanpercentile(bs, 95)), common_offset=c0,
                   med_abs_dv_station=float(ST.absdv.median()), med_err_pooled=float(np.median(np.concatenate([ST.eA, ST.eB]))),
                   by_class={k: rstat(v) for k, v in ST.groupby("dclass") if len(v) >= 5}, R_no_offset=None)
    S["F8"]["limitation_replicated"] = bool(S["F8"]["R"] < 0.5)
    dd = (ST.eA - ST.eB).values; nz = dd[dd != 0]
    S["F9"] = dict(n_B_closer=int((dd > 0).sum()), n_A_closer=int((dd < 0).sum()), sign_test_p=float(binomtest(int((nz > 0).sum()), len(nz), 0.5).pvalue) if len(nz) else None)
    S["F10"] = dict(common_offset=c0, by_class={k: float(np.median((v.vA + v.vB) / 2 - v.v_gnss)) for k, v in ST.groupby("dclass")})
else:
    S["F7"] = dict(testable=False, **{"pass": False}); S["F8"] = dict(n=nst); S["F9"] = {}; S["F10"] = {}
ST.to_csv(f"{R}/by_station.csv", index=False)
# F11
S["F11"] = {k: dict(n_tiles=int(len(v)), trend=float(v.n_trend_changed.sum() / v.n_lab.sum()), step=float(v.n_step_changed.sum() / v.n_lab.sum()),
                    m_median=float(v.m.median())) for k, v in inc.groupby("dclass")}
json.dump(S, open(f"{R}/summary.json", "w"), indent=2, ensure_ascii=False, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("written", R)
