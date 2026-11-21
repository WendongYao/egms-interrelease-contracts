"""Confirm02_Mechanism s95 (cloud): Part I (M1-M4), Part II (M5-M7), U1-U3 and the Fig. 2 rule (CONTRACT sections 3-5).
Sets: O (P1, dg == 0 by documentation), O' and N' (P2, dg = gnss_velocity_u_B - gnss_velocity_u_A per point, from s94).
Comparisons (velocities of release B harmonized, release A unchanged; both velocity variants pub and cw):
  C0 as is; C1 vB - dg; C1c vB - const (const = median over the set's included tiles of the tile-median dg); C2 = C1 then minus the
  tile least-squares plane of (vB - dg - vA); C2p = vB minus the plane of (vB - vA). Functions: rsc/rsc_core.harmonize, rsc_core.plane.
Decision metrics per tile: Confirm01 s89 tile_metrics (p2lib, verbatim slice) -> identical definitions; tile inclusion as s89.
Aggregates (as Confirm01 s89 primary): A1 = median over tiles of 'dis'; A2 = 1 - median over tiles of persA (T5, T10);
  A3 = pooled built_dis / built_flag (T5, T10).
Shares (M5): per metric, tiles with a finite, positive C0 tile value (A1: dis; A2: 1 - persA; A3: built_dis/built_flag) form the
  share subset; D_C0, D_C1, D_C2 are the aggregates over that subset; s_g = (D0-D1)/D0, s_ramp = (D1-D2)/D0, s_local = D2/D0;
  tile bootstrap 1000, numpy default_rng(20261118), 2.5/97.5 percentiles. For O, D_C1 := D_C0.
Fixed here before reading any s94 output (interpretations of the frozen contract, logged in RUN_LOG):
  - M3 is judged on common-window velocities (Confirm01 m was defined on the common-window primary model); pub reported.
    Dense = n (min epochs of the two releases) >= 200 in both O and O' for the paired O/O' comparison.
  - Fig. 2 'N' pooled value' = the N' A2 primary aggregate (median over tiles of T10 persA, pub, C0).
  - M7 B3: the B3 statistics are computed on detrended, deseasoned difference series; subtracting a constant velocity dg removes
    a linear trend and leaves them unchanged by construction, so they are reported as unchanged (not recomputed).
  - U2 floating tolerance: after C2 the recovered velocities that differ from vA by <= 1e-9 mm/yr are set to vA before statuses
    are compared (published velocities sit exactly on thresholds such as 0.5 or 5.0); the maximum residual is reported.
  - U1 also runs rsc_core.decision_metrics (generalized hotspot grid) on the pub variant for C0, C1 and C2 and requires equality
    with the p2lib values for every A1/A2/A3 field.
Usage: python s95_mech.py OUT GHSL_DIR C01_RUN_DIR G_OPRIME_DIR G_NPRIME_DIR [NPROC]"""
import sys, os, glob, json, re
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import spearmanr, rankdata
OUT, GH, RUN, G_OP, G_NP = sys.argv[1:6]; NPROC = int(sys.argv[6]) if len(sys.argv) > 6 else 2; os.makedirs(OUT, exist_ok=True)
os.environ["GHSL_DIR"] = GH
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, f"{HERE}/rsc")
import p2lib, rsc_core as rc
SEED, NBOOT = 20261118, 1000
SETS = {"O": (f"{RUN}/fulltile_O", None), "O'": (f"{RUN}/fulltile_Oprime", G_OP), "N'": (f"{RUN}/fulltile_Nprime", G_NP)}
MODES = {"O": ["C0", "C1", "C2"], "O'": ["C0", "C1", "C1c", "C2", "C2p"], "N'": ["C0", "C1", "C1c", "C2", "C2p"]}
def included(m): return ("excluded" not in m) and m["nA_epochs"] >= 60 and m["nB_epochs"] >= 60 and m["W_days"] >= 1095 and m["n_matched"] >= 50
def flat(res, pre):
    row = {}
    for (var, th), r in res.items():
        if th == "curve":
            for T, v in r.items(): row[f"{pre}{var}_curve_{T}"] = v
        else:
            for k, v in r.items(): row[f"{pre}{var}_{th}_{k}"] = v
    return row
def load_g(gdir, t):
    z = np.load(f"{gdir}/{t}.npz"); return z["gA"].astype(np.float64), z["gB"].astype(np.float64)
def dg_for(name, t, z):
    if name == "O": return np.zeros(len(z["E"]))
    gA, gB = load_g(SETS[name][1], t); return gB - gA
def dg_injection(name, t, z):
    """U2: dg_t of the tile in P2; for O taken from O' of the same tile by coordinate key (missing -> tile median)."""
    if name != "O": return dg_for(name, t, z)
    zp = np.load(SETS["O'"][0] + f"/{t}.npz"); gA, gB = load_g(G_OP, t); d = gB - gA
    kp = rc.key(zp["E"], zp["N"]); k = rc.key(z["E"], z["N"]); o = np.argsort(kp); i = np.searchsorted(kp[o], k); i = np.clip(i, 0, len(kp) - 1)
    ok = kp[o][i] == k; return np.where(ok, d[o][i], np.median(d))
def work(args):
    name, t, const = args
    fd = SETS[name][0]; m = json.load(open(f"{fd}/{t}.json"))
    if not included(m): return dict(tile=t, included=False)
    z = np.load(f"{fd}/{t}.npz"); E, N = z["E"], z["N"]; expo = p2lib.exposure(E, N); dg = dg_for(name, t, z)
    row = dict(tile=t, included=True, n=min(m["nA_epochs"], m["nB_epochs"]), Nidx=int(re.findall(r"N(\d+)", t)[0]), med_dg=float(np.median(dg)))
    vpub = z["vB_pub"].astype(float) - z["vA_pub"].astype(float); f = np.isfinite(vpub)
    row["med_dv_pub"] = float(np.median(vpub[f])); row["med_dv_minus_dg_pub"] = float(np.median((vpub - dg)[f]))
    for mode in MODES[name]:
        zz = dict(E=E, N=N)
        for var in ("pub", "cw"):
            vA, vB = z[f"vA_{var}"].astype(float), z[f"vB_{var}"].astype(float)
            vBh = rc.harmonize(E, N, vA, vB, dg, mode if mode != "C1c" else "C1c", const)
            zz[f"vA_{var}"], zz[f"vB_{var}"] = vA, vBh; d = vBh - vA; ff = np.isfinite(d)
            row[f"{mode}__mstar_{var}"] = float(np.median(np.abs(d[ff]))) if ff.any() else np.nan
            if var == "pub" and mode in ("C0", "C1", "C2"):  # U1b: RSC core on the same arrays
                r2 = rc.decision_metrics(E, N, vA, vBh, expo["built"], expo["pop"])
                for th, r in r2.items():
                    for k, v in r.items(): row[f"rsc_{mode}__pub_{th}_{k}"] = v
        row.update(flat(p2lib.tile_metrics(zz, t, expo), f"{mode}__"))
    # U2 injection
    dgi = dg_injection(name, t, z); zz = {k: dict(E=E, N=N) for k in ("I0", "I1", "I2")}
    for var in ("pub", "cw"):
        vA = z[f"vA_{var}"].astype(float); vB = z[f"vB_{var}"].astype(float); P = rc.plane(E, N, vB - dgi - vA); vBs = vA + dgi + P
        for k, mode in (("I0", "C0"), ("I1", "C1"), ("I2", "C2")):
            zz[k][f"vA_{var}"] = vA; zz[k][f"vB_{var}"] = rc.harmonize(E, N, vA, vBs, dgi, mode)
        r2 = zz["I2"][f"vB_{var}"] - vA; row[f"inj_maxabs_C2_{var}"] = float(np.nanmax(np.abs(r2)))
        # floating tolerance of U2 (contract: 1e-9) applied to the recovered velocities before statuses are compared
        zz["I2"][f"vB_{var}"] = np.where(np.abs(r2) <= 1e-9, vA, zz["I2"][f"vB_{var}"])
    for k in ("I0", "I1", "I2"): row.update(flat(p2lib.tile_metrics(zz[k], t, expo), f"{k}__"))
    return row
# ---------- run tiles
TAB = {}
for name, (fd, gd) in SETS.items():
    tiles = sorted(os.path.basename(p)[:-5] for p in glob.glob(f"{fd}/*.json"))
    const = 0.0
    if gd:
        meds = [float(np.median(np.subtract(*load_g(gd, t)[::-1]))) for t in tiles if included(json.load(open(f"{fd}/{t}.json")))]
        const = float(np.median(meds))
    with Pool(NPROC) as p: rows = p.map(work, [(name, t, const) for t in tiles], chunksize=1)
    T = pd.DataFrame(rows); T.to_csv(f"{OUT}/by_tile_{name}.csv", index=False)
    TAB[name] = (pd.read_csv(f"{OUT}/by_tile_{name}.csv"), const); print(name, "tiles", len(T), "const", const, flush=True)
# ---------- helpers
def boot_idx(n):
    g = np.random.default_rng(SEED); return [g.integers(0, n, n) for _ in range(NBOOT)]
METRICS = [("A1", th) for th in ("T0", "T05", "T5", "T10")] + [("A2", th) for th in ("T5", "T10")] + [("A3", th) for th in ("T5", "T10")]
def tileD(I, mode, var, met, th):
    p = f"{mode}__{var}_{th}_"
    if met == "A1": return I[p + "dis"].values
    if met == "A2": return 1 - I[p + "persA"].values
    return (I[p + "built_dis"] / I[p + "built_flag"]).values
def aggD(I, mode, var, met, th, idx=None):
    J = I if idx is None else I.iloc[idx]; p = f"{mode}__{var}_{th}_"
    if met == "A1": return float(np.nanmedian(J[p + "dis"]))
    if met == "A2": return float(1 - np.nanmedian(J[p + "persA"]))
    s = J[p + "built_flag"].sum(); return float(J[p + "built_dis"].sum() / s) if s > 0 else np.nan
def shares(I, var, met, th, has_c1=True):
    d0 = tileD(I, "C0", var, met, th); S = I[np.isfinite(d0) & (d0 > 0)].reset_index(drop=True)
    def f(idx):
        D0 = aggD(S, "C0", var, met, th, idx); D1 = aggD(S, "C1", var, met, th, idx) if has_c1 else D0; D2 = aggD(S, "C2", var, met, th, idx)
        return D0, D1, D2
    D0, D1, D2 = f(None); out = dict(n_tiles=int(len(S)), D0=D0, D1=D1, D2=D2, s_g=(D0 - D1) / D0, s_ramp=(D1 - D2) / D0, s_local=D2 / D0)
    bs = np.array([[(a - b) / a, (b - c) / a, c / a] for a, b, c in (f(ix) for ix in boot_idx(len(S)))])
    for j, k in enumerate(("s_g", "s_ramp", "s_local")): out[k + "_ci"] = [float(np.nanpercentile(bs[:, j], 2.5)), float(np.nanpercentile(bs[:, j], 97.5))]
    return out
def eq_frames(a, b, cols):
    bad = []
    for c in cols:
        x, y = a[c].values, b[c].values
        same = (x == y) | (pd.isna(x) & pd.isna(y))
        if not np.all(same): bad.append(c)
    return bad
R = dict(consts={k: v[1] for k, v in TAB.items()})
I = {k: v[0][v[0].included == True].reset_index(drop=True) for k, v in TAB.items()}
# ---------- U1 / U3
c01 = {"O": "by_tile_O.csv", "O'": "by_tile_O'.csv", "N'": "by_tile_N'.csv"}; U1 = {}
for name in SETS:
    ref = pd.read_csv(f"{RUN}/part2/{c01[name]}"); mine = TAB[name][0]
    cols = [c for c in ref.columns if c not in ("tile", "included", "n", "Nidx")]
    a = mine.set_index("tile").rename(columns={f"C0__{c}": c for c in cols}); b = ref.set_index("tile")
    common = [c for c in cols if c in a.columns]; b = b.loc[[t for t in b.index if t in a.index]]; a = a.loc[b.index]
    U1[name] = dict(n_tiles_ref=int(len(ref)), n_tiles_compared=int(len(b)), n_cols=len(common), n_cols_missing=len(cols) - len(common),
                    mismatched_cols=eq_frames(a.reset_index(), b.reset_index(), common))
    rb = []
    for mode in ("C0", "C1", "C2"):
        for th in ("T5", "T10", "T0", "T05"):
            for k in ("dis", "n_dis", "n_flag", "built_dis", "built_flag", "pop_dis", "pop_flag") + (("loss", "gain", "n_hotA", "n_hotB", "persA", "persB", "iou_med") if th in ("T5", "T10") else ()):
                c1, c2 = f"rsc_{mode}__pub_{th}_{k}", f"{mode}__pub_{th}_{k}"
                if c1 in mine and c2 in mine:
                    x, y = mine[c1].values, mine[c2].values
                    if not np.all((x == y) | (pd.isna(x) & pd.isna(y))): rb.append(c1)
    U1[name]["rsc_vs_p2lib_mismatched"] = rb
U1["pass"] = bool(all(not U1[n]["mismatched_cols"] and not U1[n]["rsc_vs_p2lib_mismatched"] and U1[n]["n_cols_missing"] == 0 and U1[n]["n_tiles_ref"] == U1[n]["n_tiles_compared"] for n in SETS))
R["U1"] = U1
c1cols = [c for c in TAB["O"][0].columns if c.startswith("C1__")]
R["U3"] = dict(mismatched=eq_frames(TAB["O"][0].rename(columns=lambda c: c.replace("C1__", "X__")), TAB["O"][0].rename(columns=lambda c: c.replace("C0__", "X__")), [c.replace("C1__", "X__") for c in c1cols]))
R["U3"]["pass"] = not R["U3"]["mismatched"]
# ---------- U2
U2 = {}
for name in SETS:
    J = I[name]; worst = 0.0; s_ok = True
    for var in ("pub", "cw"):
        for met, th in METRICS:
            v2 = tileD(J, "I2", var, met, th) if met != "A2" else (1 - J[f"I2__{var}_{th}_persA"].values)
            if met == "A3": v2 = J[f"I2__{var}_{th}_built_dis"].values
            worst = max(worst, float(np.nanmax(np.abs(np.nan_to_num(v2)))))
        # whole-set aggregate share identity on A1 (all thresholds)
        for th in ("T0", "T05", "T5", "T10"):
            D0 = float(np.nanmedian(J[f"I0__{var}_{th}_dis"])); D1 = float(np.nanmedian(J[f"I1__{var}_{th}_dis"])); D2 = float(np.nanmedian(J[f"I2__{var}_{th}_dis"]))
            if D0 > 0 and not (abs((D0 - D1) / D0 + (D1 - D2) / D0 - 1) <= 1e-9 and abs(D2 / D0) <= 1e-9): s_ok = False
    U2[name] = dict(max_abs_D_after_C2=worst, max_abs_velocity_residual=float(J[["inj_maxabs_C2_pub", "inj_maxabs_C2_cw"]].max().max()), share_identity_ok=s_ok)
U2["pass"] = bool(all(U2[n]["max_abs_D_after_C2"] <= 1e-9 and U2[n]["share_identity_ok"] for n in SETS)); R["U2"] = U2
# ---------- M1, M2
for name in ("O'", "N'"):
    J = I[name]; md = J.med_dg.values
    R.setdefault("M1", {})[name] = dict(n=int(len(md)), frac_in_range=float(np.mean((md >= -0.11) & (md <= 0.98))), frac_pos=float(np.mean(md > 0)),
                                         median=float(np.median(md)), p05=float(np.percentile(md, 5)), p95=float(np.percentile(md, 95)), min=float(md.min()), max=float(md.max()))
    R["M1"][name]["holds"] = bool(R["M1"][name]["frac_in_range"] >= 0.95 and R["M1"][name]["frac_pos"] >= 0.90)
    S0 = float(np.median(J.med_dv_pub)); S1 = float(np.median(J.med_dv_minus_dg_pub))
    rho = float(spearmanr(J.med_dv_pub, J.med_dg).statistic); sl = float(np.polyfit(J.med_dg, J.med_dv_pub, 1)[0]) if J.med_dg.std() > 0 else np.nan
    bs_r, bs_s = [], []
    for ix in boot_idx(len(J)):
        s = J.iloc[ix]; bs_r.append(spearmanr(s.med_dv_pub, s.med_dg).statistic)
        bs_s.append(np.polyfit(s.med_dg, s.med_dv_pub, 1)[0] if s.med_dg.std() > 0 else np.nan)
    R.setdefault("M2", {})[name] = dict(S0=S0, S1=S1, holds=bool(abs(S1) <= 0.1 and abs(S1) <= abs(S0) / 3), rho=rho,
                                         rho_ci=[float(np.nanpercentile(bs_r, 2.5)), float(np.nanpercentile(bs_r, 97.5))], slope=sl,
                                         slope_ci=[float(np.nanpercentile(bs_s, 2.5)), float(np.nanpercentile(bs_s, 97.5))])
R["M1"]["holds_both"] = bool(R["M1"]["O'"]["holds"] and R["M1"]["N'"]["holds"]); R["M2"]["holds_both"] = bool(R["M2"]["O'"]["holds"] and R["M2"]["N'"]["holds"])
# ---------- M3
O_, Op = I["O"].set_index("tile"), I["O'"].set_index("tile"); both = sorted(set(O_[O_.n >= 200].index) & set(Op[Op.n >= 200].index))
M3 = dict(n_paired_dense=len(both))
for var in ("cw", "pub"):
    mo = float(np.median(O_.loc[both, f"C0__mstar_{var}"])); m0 = float(np.median(Op.loc[both, f"C0__mstar_{var}"])); m1 = float(np.median(Op.loc[both, f"C1__mstar_{var}"]))
    M3[var] = dict(O_C0=mo, Oprime_C0=m0, Oprime_C1=m1, ratio_C0=m0 / mo, ratio_C1=m1 / mo)
for name in SETS:
    J = I[name][I[name].n >= 200]
    M3[f"all_dense_{name}"] = {f"{mode}_{var}": float(np.median(J[f"{mode}__mstar_{var}"])) for mode in MODES[name] for var in ("cw", "pub")}
M3["holds"] = bool(2 / 3 <= M3["cw"]["ratio_C1"] <= 1.5); R["M3"] = M3
# ---------- M5, M6
M5 = {}
for name in SETS:
    for var in ("pub", "cw"):
        for met, th in METRICS: M5[f"{name}|{var}|{met}|{th}"] = shares(I[name], var, met, th, has_c1=(name != "O"))
a_ok = all(M5[f"{n}|pub|A1|{th}"]["s_g"] >= 0.3 and M5[f"{n}|pub|A1|{th}"]["s_g_ci"][0] > 0 for n in ("O'", "N'") for th in ("T0", "T05"))
b_ok = all(M5[f"{n}|pub|{met}|T10"]["s_g_ci"][1] <= 0.10 for n in ("O'", "N'") for met in ("A1", "A2", "A3"))
cO = aggD(I["O"], "C0", "pub", "A1", "T05"); cP = aggD(I["O'"], "C1", "pub", "A1", "T05"); cP0 = aggD(I["O'"], "C0", "pub", "A1", "T05")
R["M5"] = dict(shares=M5, a_holds=bool(a_ok), b_holds=bool(b_ok), c=dict(O_C0=cO, Oprime_C0=cP0, Oprime_C1=cP, ratio_C1=cP / cO, holds=bool(2 / 3 <= cP / cO <= 1.5)))
M6 = {}
for name in ("O'", "N'"):
    for var in ("pub", "cw"):
        for met, th in METRICS:
            M6[f"{name}|{var}|{met}|{th}"] = dict(C1=aggD(I[name], "C1", var, met, th), C1c=aggD(I[name], "C1c", var, met, th),
                                                 C2=aggD(I[name], "C2", var, met, th), C2p=aggD(I[name], "C2p", var, met, th))
            M6[f"{name}|{var}|{met}|{th}"]["abs_C1c_minus_C1"] = abs(M6[f"{name}|{var}|{met}|{th}"]["C1c"] - M6[f"{name}|{var}|{met}|{th}"]["C1"])
            M6[f"{name}|{var}|{met}|{th}"]["abs_C2p_minus_C2"] = abs(M6[f"{name}|{var}|{met}|{th}"]["C2p"] - M6[f"{name}|{var}|{met}|{th}"]["C2"])
R["M6"] = M6
# also: full-set aggregates of every comparison (for tables)
R["aggregates"] = {f"{n}|{mode}|{var}|{met}|{th}": aggD(I[n], mode, var, met, th) for n in SETS for mode in MODES[n] for var in ("pub", "cw") for met, th in METRICS}
# ---------- M4 (GNSS ratio R under C1) and M7 (B1 under C1)
def rstat(d):
    c = float(np.median((d.vA + d.vB) / 2 - d.v_gnss)); eA = np.abs(d.vA - c - d.v_gnss).values; eB = np.abs(d.vB - c - d.v_gnss).values
    den = float(np.median(np.concatenate([eA, eB]))); return float(np.median(np.abs(d.vB - d.vA))) / den if den > 0 else np.nan
def pcorr(x, y, zc):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(zc); X = np.vstack([np.ones_like(rz), rz]).T
    ex = rx - X @ np.linalg.lstsq(X, rx, rcond=None)[0]; ey = ry - X @ np.linalg.lstsq(X, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])
def boot90(df, f):
    g = np.random.default_rng(20261117); out = []
    for _ in range(1000):
        s = df.iloc[g.integers(0, len(df), len(df))]
        try: out.append(f(s))
        except Exception: out.append(np.nan)
    return float(np.nanpercentile(out, 2.5)), float(np.nanpercentile(out, 97.5))
M4, M7 = {}, {}
for name, sd, gd in (("O'", "Oprime", G_OP), ("N'", "Nprime", G_NP)):
    P = pd.read_csv(f"{RUN}/part1_{sd}/by_point.csv"); fd = SETS[name][0]; P["dg"] = np.nan
    P = P[[os.path.exists(f"{gd}/{t}.npz") for t in P.tile]].copy()
    for t, idx in P.groupby("tile").groups.items():
        z = np.load(f"{fd}/{t}.npz"); gA, gB = load_g(gd, t); kf = rc.key(z["E"], z["N"]); o = np.argsort(kf)
        kp = rc.key(P.loc[idx, "E"].values, P.loc[idx, "N"].values); i = np.clip(np.searchsorted(kf[o], kp), 0, len(kf) - 1); ok = kf[o][i] == kp
        P.loc[idx, "dg"] = np.where(ok, (gB - gA)[o][i], np.nan)
    M4[name] = dict(n_points_sta=int((P.kind == "sta").sum()), n_points_sta_unmatched=int(((P.kind == "sta") & P.dg.isna()).sum()))
    S = P[(P.kind == "sta") & P.dg.notna()]; G = pd.read_csv(f"{RUN}/gnss_{sd}/gnss_fit.csv"); G = G[G.included == True][["sta", "v_gnss"]]
    st = []
    for (s, t), d in S.groupby(["sta", "tile"]):
        if len(d) >= 3: st.append(dict(sta=s, tile=t, vA=float(np.median(d.vA)), vB=float(np.median(d.vB)), vB1=float(np.median(d.vB - d.dg))))
    ST = pd.DataFrame(st, columns=["sta", "tile", "vA", "vB", "vB1"]).merge(G, on="sta", how="inner"); ST1 = ST.assign(vB=ST.vB1)
    if len(ST) < 5:
        M4[name]["untestable"] = True; M7[name] = dict(untestable=True); continue
    g = np.random.default_rng(SEED); bsi = [g.integers(0, len(ST), len(ST)) for _ in range(NBOOT)]
    M4[name].update(n_stations=int(len(ST)), R_C0=rstat(ST), R_C1=rstat(ST1),
                    R_C1_p05_p95=[float(np.nanpercentile([rstat(ST1.iloc[ix]) for ix in bsi], 5)), float(np.nanpercentile([rstat(ST1.iloc[ix]) for ix in bsi], 95))])
    Rr = P[P.kind == "rand"]; T = pd.read_csv(f"{RUN}/part1_{sd}/by_tile.csv"); T = T[T.included == True].copy()
    mC0 = Rr.groupby("tile").apply(lambda d: float(np.median(np.abs(d.vB - d.vA)))); mC1 = Rr[Rr.dg.notna()].groupby("tile").apply(lambda d: float(np.median(np.abs(d.vB - d.dg - d.vA))))
    T = T[T.tile.isin(mC1.index)].copy(); T["m_C0"] = T.tile.map(mC0); T["m_C1"] = T.tile.map(mC1); T["G"] = [max(json.load(open(f"{fd}/{t}.json"))["gapA_days"], json.load(open(f"{fd}/{t}.json"))["gapB_days"]) for t in T.tile]
    T = T[T.m_C1.notna()]
    if len(T) < 5:
        M7[name] = dict(untestable=True); continue
    r1, r2 = pcorr(T.m_C1, T.G, T.n), pcorr(T.m_C1, T.n, T.G); ci1, ci2 = boot90(T, lambda s: pcorr(s.m_C1, s.G, s.n)), boot90(T, lambda s: pcorr(s.m_C1, s.n, s.G))
    ex0 = lambda ci: ci[0] <= 0 <= ci[1]
    v = "winter gap dominates" if (abs(r1) >= 0.3 and not ex0(ci1) and ex0(ci2)) else ("sampling dominates" if (abs(r2) >= 0.3 and not ex0(ci2) and ex0(ci1)) else "undetermined")
    M7[name] = dict(B1_C1=dict(rho_m_G_given_n=r1, ci_m_G_given_n=ci1, rho_m_n_given_G=r2, ci_m_n_given_G=ci2, verdict=v),
                    m_C0_equals_part1_m=bool(np.allclose(T.m_C0.values, T.m.values, rtol=0, atol=1e-12)),
                    m_dense_C0=float(T[T.dclass == "dense"].m_C0.median()), m_dense_C1=float(T[T.dclass == "dense"].m_C1.median()),
                    m_sparse_C0=float(T[T.dclass == "sparse"].m_C0.median()) if (T.dclass == "sparse").any() else None,
                    m_sparse_C1=float(T[T.dclass == "sparse"].m_C1.median()) if (T.dclass == "sparse").any() else None,
                    B3="unchanged by construction (detrended difference series; dg is a constant velocity)")
R["M4"], R["M7"] = M4, M7
# ---------- Fig. 2 rule
J = I["N'"]; agg = float(np.nanmedian(J["C0__pub_T10_persA"])); C = J[J["C0__pub_T10_n_hotA"] >= 10].copy()
C["dist"] = (C["C0__pub_T10_persA"] - agg).abs(); C = C.sort_values(["dist", "tile"])
R["Fig2"] = dict(pooled=agg, n_candidates=int(len(C)), tile=C.tile.iloc[0] if len(C) else None, persA=float(C["C0__pub_T10_persA"].iloc[0]) if len(C) else None)
json.dump(R, open(f"{OUT}/summary.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("written", OUT)
