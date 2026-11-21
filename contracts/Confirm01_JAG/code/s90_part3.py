"""Confirm01_JAG s90 (cloud): Part III (B) attribution, CONTRACT section 4 Part III and Amendment 1 section 2.4.
Primary set N' (new tiles x P2); control O' (old tiles x P2). Inputs: Part I outputs of s88 (by_tile.csv: m, n, dclass),
full-tile outputs of s82 (meta: longest gaps, epochs; npz: height, ac1_diff, rescorr) and s86 samples (for B2), and the
N' manifest (stratum, for the B2 donor pool).
Tile covariates: G = max(gapA_days, gapB_days) inside W; n; mean northing of the sample points; median point height.
B1: partial Spearman rho(m, G | n) and rho(m, n | G) (ranks; linear regression of ranks on the conditioning rank; correlation
    of residuals); bootstrap over tiles, 1000, seed 20261117, 95% percentile interval.
    'winter gap dominates' if |rho(m,G|n)| >= 0.3 with CI excluding 0 and the CI of rho(m,n|G) including 0; 'sampling dominates'
    if the reverse; otherwise 'undetermined'.
B2: each dense tile (n >= 200) of the set gets one donor drawn with numpy default_rng([20261117, 2]) from the set's tiles that
    are sparse (n <= 130) or in stratum Al (dense tiles excluded from the pool), in sorted tile order; both releases keep only
    epochs within 3 days of one of the donor's own epochs of the same release inside the donor's W; m is recomputed on the same
    250 sample points (primary model, min 60 epochs). S = log(median over dense tiles of m_thin/m_full) / log(median m of sparse
    tiles / median m of dense tiles); bootstrap over dense tiles for the interval. S <= 0.25 with upper < 0.5: sampling explains a
    minor part; S >= 0.5: sampling is the main factor; otherwise undetermined.
B3: median over tiles of the tile-median lag-1 autocorrelation of the detrended, deseasoned B-A difference series >= 0.4 ->
    'low-frequency ramp' holds; median over dense tiles minus median over sparse tiles of the tile-median residual correlation >= 0.15
    -> holds.
Usage: python s90_part3.py OUT_DIR SETNAME PART1_DIR FULLTILE_DIR SAMPLE_DIR MANIFEST [SETNAME PART1_DIR FULLTILE_DIR SAMPLE_DIR MANIFEST ...]"""
import sys, os, glob, json, math
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rc_common import window, tyears, fit_primary, d2ord
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True); A = sys.argv[2:]; SEED = 20261117; NBOOT = 1000
def pcorr(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z); X = np.vstack([np.ones_like(rz), rz]).T
    ex = rx - X @ np.linalg.lstsq(X, rx, rcond=None)[0]; ey = ry - X @ np.linalg.lstsq(X, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])
def boot(df, f):
    g = np.random.default_rng(SEED); out = []
    for _ in range(NBOOT):
        s = df.iloc[g.integers(0, len(df), len(df))]
        try: out.append(f(s))
        except Exception: out.append(np.nan)
    return float(np.nanpercentile(out, 2.5)), float(np.nanpercentile(out, 97.5))
summ = {}
for k in range(0, len(A), 5):
    name, P1, FD, SD, MAN = A[k:k + 5]
    T = pd.read_csv(f"{P1}/by_tile.csv"); T = T[T.included == True].copy(); man = pd.read_csv(MAN).set_index("tile")
    cov = []
    for t in T.tile:
        m = json.load(open(f"{FD}/{t}.json")); z = np.load(f"{FD}/{t}.npz")
        cov.append(dict(tile=t, G=max(m["gapA_days"], m["gapB_days"]), height_med=float(np.nanmedian(z["height"])) if np.isfinite(z["height"]).any() else np.nan,
                        ac1_med=float(np.nanmedian(z["ac1_diff"])), rescorr_med=float(np.nanmedian(z["rescorr"]))))
    T = T.merge(pd.DataFrame(cov), on="tile"); T["stratum"] = [man.loc[t, "stratum"] if "stratum" in man.columns else "" for t in T.tile]
    S = dict(n_tiles=int(len(T)))
    # B1
    r1, r2 = pcorr(T.m, T.G, T.n), pcorr(T.m, T.n, T.G)
    ci1, ci2 = boot(T, lambda s: pcorr(s.m, s.G, s.n)), boot(T, lambda s: pcorr(s.m, s.n, s.G))
    ex0 = lambda ci: ci[0] <= 0 <= ci[1]
    verdict = "winter gap dominates" if (abs(r1) >= 0.3 and not ex0(ci1) and ex0(ci2)) else ("sampling dominates" if (abs(r2) >= 0.3 and not ex0(ci2) and ex0(ci1)) else "undetermined")
    S["B1"] = dict(rho_m_G_given_n=r1, ci_m_G_given_n=ci1, rho_m_n_given_G=r2, ci_m_n_given_G=ci2, verdict=verdict,
                   rho_m_northing=float(spearmanr(T.m, T.northing_mean).statistic), rho_m_height=float(spearmanr(T.m, T.height_med, nan_policy="omit").statistic),
                   rho_G_n=float(spearmanr(T.G, T.n).statistic))
    # B2
    dense = T[T.dclass == "dense"].sort_values("tile"); pool = T[(T.dclass == "sparse") | (T.stratum == "Al")]; pool = pool[pool.dclass != "dense"].sort_values("tile")
    rows = []
    if len(dense) and len(pool):
        g = np.random.default_rng([SEED, 2]); donors = pool.tile.values[g.integers(0, len(pool), len(dense))]
        for t, dn in zip(dense.tile, donors):
            zt = np.load(f"{SD}/{t}.npz", allow_pickle=True); zd = np.load(f"{SD}/{dn}.npz", allow_pickle=True)
            lo, hi, ia, ib = window(zt["dA"], zt["dB"]); dA = [str(zt["dA"][i]) for i in ia]; dB = [str(zt["dB"][i]) for i in ib]
            lo2, hi2, ja, jb = window(zd["dA"], zd["dB"]); oA = np.array([d2ord(str(zd["dA"][i])) for i in ja]); oB = np.array([d2ord(str(zd["dB"][i])) for i in jb])
            keepA = [i for i, d in enumerate(dA) if np.min(np.abs(oA - d2ord(d))) <= 3]; keepB = [i for i, d in enumerate(dB) if np.min(np.abs(oB - d2ord(d))) <= 3]
            Af, Bf = zt["A_rand"][:, ia].astype(float), zt["B_rand"][:, ib].astype(float); tA, tB = tyears(dA), tyears(dB)
            full, thin = [], []
            for p in range(Af.shape[0]):
                a, b = fit_primary(tA, Af[p]), fit_primary(tB, Bf[p])
                a2, b2 = fit_primary(tA[keepA], Af[p, keepA]), fit_primary(tB[keepB], Bf[p, keepB])
                if a and b: full.append(abs(b[0] - a[0]))
                if a2 and b2: thin.append(abs(b2[0] - a2[0]))
            rows.append(dict(tile=t, donor=dn, n_keepA=len(keepA), n_keepB=len(keepB), m_full=float(np.median(full)) if full else np.nan,
                             m_thin=float(np.median(thin)) if len(thin) >= 50 else np.nan, n_thin_points=len(thin)))
    B2 = pd.DataFrame(rows); B2.to_csv(f"{OUT}/B2_{name}.csv", index=False)
    md, ms = T[T.dclass == "dense"].m.median(), T[T.dclass == "sparse"].m.median()
    if len(B2) and B2.m_thin.notna().sum() >= 5 and ms > md > 0:
        def Sfun(d): return math.log(np.nanmedian(d.m_thin / d.m_full)) / math.log(ms / md)
        Sv = Sfun(B2); ci = boot(B2[B2.m_thin.notna()], Sfun)
        v = "minor" if (Sv <= 0.25 and ci[1] < 0.5) else ("main" if Sv >= 0.5 else "undetermined")
        S["B2"] = dict(S=Sv, ci=ci, n_dense=int(B2.m_thin.notna().sum()), median_ratio=float(np.nanmedian(B2.m_thin / B2.m_full)), m_dense=float(md), m_sparse=float(ms), verdict=v)
    else:
        S["B2"] = dict(untestable=True, n_dense=int(len(B2)), m_dense=float(md) if pd.notna(md) else None, m_sparse=float(ms) if pd.notna(ms) else None)
    # B3
    a1 = float(T.ac1_med.median()); dres = float(T[T.dclass == "dense"].rescorr_med.median() - T[T.dclass == "sparse"].rescorr_med.median()) if ((T.dclass == "dense").any() and (T.dclass == "sparse").any()) else np.nan
    S["B3"] = dict(ac1_median=a1, ramp_holds=bool(a1 >= 0.4), rescorr_by_class={c: float(v.rescorr_med.median()) for c, v in T.groupby("dclass")},
                   rescorr_dense_minus_sparse=dres, rescorr_holds=bool(dres >= 0.15) if np.isfinite(dres) else None)
    T.to_csv(f"{OUT}/by_tile_{name}.csv", index=False); summ[name] = S
json.dump(summ, open(f"{OUT}/summary.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("written", OUT)
