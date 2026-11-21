"""Confirm03_Proximity s98 (cloud): part D of CONTRACT_CN.md (SHA 257ab6ea...8af1), D-1 east-west component and D-2 sparse supplement S''.
Inputs (extracted on the user's VM with unchanged code; transfer hashes in RUN_LOG):
  D-1: Confirm01 s82 run on the E files of O' and N' (fulltile_E_*), Confirm03 s97 gnss_velocity_e (gE_*); dg_E = gE_B - gE_A.
  D-2: N' (U) = Confirm01 fulltile_Nprime + Confirm02 g_Nprime; S'' (U) = s82 on manifest_Spp + Confirm02 s94 (g_Spp).
Definitions are reused, not rewritten: decision metrics = Confirm01 p2lib.tile_metrics (via Confirm02 p2lib.py); harmonization =
  rsc_core.harmonize (C0, C1 = vB - dg, C2 = C1 minus tile plane); A1/A2/A3 aggregates and s_g/s_ramp/s_local shares = the functions
  flat, boot_idx, tileD, aggD, shares taken verbatim from Confirm02 s95_mech.py (exec of their source, SEED 20261118); part-A
  quantities (bins, N1, N2, rho1, rho2, pooled loss) = partA / summarize / ratio_ci imported from Confirm03 s96 (seed 20261120);
  M2 S0/S1, Spearman rho and slope with bootstrap = the s95 M2 block; F2 = the Confirm01 s88 F2 block (seed 20261117).
Tile inclusion as s95/s96: not excluded, both releases >= 60 epochs in W, W >= 1095 days, >= 50 matched points.
Interpretations fixed here, before any part-D output was read (also logged in RUN_LOG):
  1. H8 tile m: Confirm01 F2 used the median |vB - vA| of the random sample points (common-window primary model). S'' has no
     sample extraction, so for N' and S'' alike m = median over all full-tile points of |vB_cw - vA_cw| (C1: vB_cw - dg), the
     same full-tile quantity Confirm02 M3 used. Density class from n = min(nA_epochs, nB_epochs) with rc_common.density_class
     (dense >= 200, sparse <= 130). For N' the Spearman correlation between this m and the Confirm01 sample m is reported.
     Judgement on C1 cw; C0 cw and both pub variants reported.
  2. H7: evaluated per set. If |S0_E| >= 0.2 in both sets, H7 holds iff both satisfy |S1_E| <= 0.1 and |S1_E| <= |S0_E|/3; if in
     neither set, "untestable"; if in one set only, that set alone is judged and the result is labelled "one_set_only".
  3. "T5/T10 loss and gain" is reported both as the median over tiles of the p2lib tile loss/gain and as the pooled part-A rate.
  4. N' sparse class = included N' tiles with n <= 130.
  5. M1_E is reported only (the U range criterion of Confirm02 M1 is not applied to E). T0 is not reported for E (contract D-1).
Outputs: OUT/summary_D.json, OUT/D_tile_E.csv, OUT/D_tile_U.csv, OUT/D_A_tile.csv. Prints progress and counts only.
Usage: python s98_confirm03_D.py OUT DATAD C01_RUN C02_DIR GHSL_DIR [NPROC] [SMOKE tiles comma-separated]"""
import sys, os, glob, json, re, ast
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import spearmanr
OUT, DD, RUN, C02, GH = sys.argv[1:6]; NPROC = int(sys.argv[6]) if len(sys.argv) > 6 else 2
SMOKE = sys.argv[7].split(",") if len(sys.argv) > 7 else None
os.makedirs(OUT, exist_ok=True); os.environ["GHSL_DIR"] = GH
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, f"{C02}/code"); sys.path.insert(0, f"{C02}/code/rsc"); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(RUN.rstrip("/")), "code"))
import p2lib, rsc_core as rc
from rc_common import density_class
_argv = sys.argv; sys.argv = ["s96", os.path.join(OUT, "_s96_unused"), RUN, C02, DD]; import s96_confirm03 as s96; sys.argv = _argv
# verbatim functions from Confirm02 s95
SEED, NBOOT = 20261118, 1000
_src = open(f"{C02}/code/s95_mech.py").read(); _tree = ast.parse(_src)
for node in _tree.body:
    if isinstance(node, ast.FunctionDef) and node.name in ("flat", "boot_idx", "tileD", "aggD", "shares"):
        exec(compile(ast.Module([node], []), "s95_mech.py", "exec"), globals())

SETS = {"O'E": (f"{DD}/fulltile_E_Oprime", f"{DD}/gE_Oprime"), "N'E": (f"{DD}/fulltile_E_Nprime", f"{DD}/gE_Nprime"),
        "N'": (f"{RUN}/fulltile_Nprime", f"{C02}/data/g_Nprime"), "S''": (f"{DD}/fulltile_Spp", f"{DD}/g_Spp")}
MODES = ["C0", "C1", "C2"]
def included(m): return ("excluded" not in m) and m["nA_epochs"] >= 60 and m["nB_epochs"] >= 60 and m["W_days"] >= 1095 and m["n_matched"] >= 50

def work(args):
    name, t = args; fd, gd = SETS[name]; m = json.load(open(f"{fd}/{t}.json"))
    if not included(m): return dict(set=name, tile=t, included=False), []
    z = np.load(f"{fd}/{t}.npz"); E, N = z["E"], z["N"]; g = np.load(f"{gd}/{t}.npz")
    dg = g["gB"].astype(np.float64) - g["gA"].astype(np.float64); expo = p2lib.exposure(E, N)
    n = min(m["nA_epochs"], m["nB_epochs"])
    row = dict(set=name, tile=t, included=True, n=n, dclass=density_class(n), med_dg=float(np.median(dg)))
    vpub = z["vB_pub"].astype(float) - z["vA_pub"].astype(float); f = np.isfinite(vpub)
    row["med_dv_pub"] = float(np.median(vpub[f])); row["med_dv_minus_dg_pub"] = float(np.median((vpub - dg)[f]))
    arows = []
    for mode in MODES:
        zz = dict(E=E, N=N)
        for var in ("pub", "cw"):
            vA, vB = z[f"vA_{var}"].astype(float), z[f"vB_{var}"].astype(float)
            vBh = rc.harmonize(E, N, vA, vB, dg, mode); zz[f"vA_{var}"], zz[f"vB_{var}"] = vA, vBh; d = vBh - vA; ff = np.isfinite(d)
            row[f"{mode}__mstar_{var}"] = float(np.median(np.abs(d[ff]))) if ff.any() else np.nan
            if mode in ("C0", "C1"):
                for T in (5.0, 10.0):
                    for dirn, ref, oth in (("loss", vA, vBh), ("gain", vBh, vA)):
                        r = s96.partA(ref, oth, T); r.update(set=name, tile=t, var=var, mode=mode, T=int(T), dir=dirn); arows.append(r)
        row.update(flat(p2lib.tile_metrics(zz, t, expo), f"{mode}__"))
    return row, arows

def m1(J):
    md = J.med_dg.values
    return dict(n=int(len(md)), median=float(np.median(md)), p05=float(np.percentile(md, 5)), p95=float(np.percentile(md, 95)),
                min=float(md.min()), max=float(md.max()), frac_pos=float(np.mean(md > 0)), frac_neg=float(np.mean(md < 0)))
def m2(J):
    S0 = float(np.median(J.med_dv_pub)); S1 = float(np.median(J.med_dv_minus_dg_pub))
    rho = float(spearmanr(J.med_dv_pub, J.med_dg).statistic); sl = float(np.polyfit(J.med_dg, J.med_dv_pub, 1)[0]) if J.med_dg.std() > 0 else np.nan
    bs_r, bs_s = [], []
    for ix in boot_idx(len(J)):
        s = J.iloc[ix]; bs_r.append(spearmanr(s.med_dv_pub, s.med_dg).statistic)
        bs_s.append(np.polyfit(s.med_dg, s.med_dv_pub, 1)[0] if s.med_dg.std() > 0 else np.nan)
    return dict(S0=S0, S1=S1, rho=rho, rho_ci=[float(np.nanpercentile(bs_r, 2.5)), float(np.nanpercentile(bs_r, 97.5))], slope=sl,
                slope_ci=[float(np.nanpercentile(bs_s, 2.5)), float(np.nanpercentile(bs_s, 97.5))])
def f2(inc, col):
    """Confirm01 s88 F2 block on a tile table with columns n, dclass and the m column `col`."""
    inc = inc[np.isfinite(inc[col])].reset_index(drop=True)
    rho = spearmanr(inc.n, inc[col]).statistic; g = np.random.default_rng(20261117); bs = []
    for _ in range(1000):
        s = inc.iloc[g.integers(0, len(inc), len(inc))]; bs.append(spearmanr(s.n, s[col]).statistic if s.n.nunique() > 1 else np.nan)
    md, ms = inc[inc.dclass == "dense"][col], inc[inc.dclass == "sparse"][col]; ratio_ok = len(md) >= 10 and len(ms) >= 10
    o = dict(n_tiles=int(len(inc)), n_by_class={k: int(v) for k, v in inc.dclass.value_counts().items()}, rho=float(rho),
             rho_ci_lo=float(np.nanpercentile(bs, 2.5)), rho_ci_hi=float(np.nanpercentile(bs, 97.5)),
             med_by_class={k: float(v.median()) for k, v in inc.groupby("dclass")[col]},
             ratio=float(ms.median() / md.median()) if (len(md) and len(ms)) else None, ratio_evaluated=bool(ratio_ok))
    o["pass"] = bool(rho <= -0.5 and o["rho_ci_hi"] < -0.2 and ((o["ratio"] >= 3) if ratio_ok else True))
    return o
def h123(AG, n, mode):
    k = lambda d: f"{n}|pub|{mode}|T5|{d}"; L, G = AG[k("loss")], AG[k("gain")]
    lab = lambda o: "above" if o["rho1_ci"][0] > 1 else ("below" if o["rho1_ci"][1] < 1 else "includes_1")
    both_in = all(0.8 <= o["rho1"] <= 1.25 and o["rho1_ci"][0] <= 1 <= o["rho1_ci"][1] for o in (L, G))
    v = "both_above" if (lab(L) == "above" and lab(G) == "above") else ("consistent_with_random_sign" if both_in else ("directional" if (lab(L) == "above") != (lab(G) == "above") else "other"))
    return dict(H1=dict(share=L["share_lost_d_lt_1"], ratio=L["rate_b0_over_rate_b3"], meets=bool(L["share_lost_d_lt_1"] >= 0.5 and L["rate_b0_over_rate_b3"] >= 3)),
                H2=dict(rho1_loss=L["rho1"], ci_loss=L["rho1_ci"], rho1_gain=G["rho1"], ci_gain=G["rho1_ci"], verdict=v),
                H3=dict(rho2=L["rho2"], ci=L["rho2_ci"], meets=bool(L["rho2"] >= 1.5 and L["rho2_ci"][0] > 1)))
def lossgain(J, A, n, mode, var, T):
    p = f"{mode}__{var}_T{T}_"
    pool = lambda d: (lambda s: float(s.n_lost.sum() / s.n_mov.sum()) if s.n_mov.sum() else np.nan)(A[(A.set == n) & (A["mode"] == mode) & (A["var"] == var) & (A["T"] == T) & (A.dir == d) & A.tile.isin(J.tile)])
    return dict(loss_tile_median=float(np.nanmedian(J[p + "loss"])), gain_tile_median=float(np.nanmedian(J[p + "gain"])), loss_pooled=pool("loss"), gain_pooled=pool("gain"))
def block(J, A, n, ths_A1, has_c1=True):
    """Report-only metric block for one tile subset J (included tiles) of set n."""
    o = dict(n_tiles=int(len(J)), aggregates={}, shares={}, lossgain={})
    for mode in MODES:
        for var in ("pub", "cw"):
            for met, th in [("A1", t) for t in ths_A1] + [("A2", "T5"), ("A2", "T10"), ("A3", "T5"), ("A3", "T10")]:
                o["aggregates"][f"{mode}|{var}|{met}|{th}"] = aggD(J, mode, var, met, th)
    for var in ("pub", "cw"):
        for met, th in [("A1", t) for t in ths_A1] + [("A2", "T5"), ("A2", "T10"), ("A3", "T5"), ("A3", "T10")]:
            try: o["shares"][f"{var}|{met}|{th}"] = shares(J.reset_index(drop=True), var, met, th, has_c1=has_c1)
            except Exception as e: o["shares"][f"{var}|{met}|{th}"] = dict(error=type(e).__name__)
        for mode in ("C0", "C1"):
            for T in (5, 10): o["lossgain"][f"{mode}|{var}|T{T}"] = lossgain(J, A, n, mode, var, T)
    AG = {}
    sub = A[A.tile.isin(J.tile) & (A.set == n)]
    for (var, mode, T, dirn), g in sub.groupby(["var", "mode", "T", "dir"]):
        r = s96.summarize(g); r.pop("bins", None); AG[f"{n}|{var}|{mode}|T{T}|{dirn}"] = r
    o["partA"] = AG; o["H123_C1"] = h123(AG, n, "C1"); o["H123_C0"] = h123(AG, n, "C0")
    return o

if __name__ == "__main__":
    rows, arows = [], []
    for name, (fd, gd) in SETS.items():
        tiles = sorted(os.path.basename(p)[:-5] for p in glob.glob(f"{fd}/*.json"))
        if SMOKE: tiles = [t for t in tiles if t in SMOKE] or tiles[:2]
        with Pool(NPROC) as p: res = p.map(work, [(name, t) for t in tiles], chunksize=1)
        for r, a in res: rows.append(r); arows += a
        print(name, "tiles", len(tiles), "included", sum(r["included"] for r, _ in res), flush=True)
    T = pd.DataFrame(rows); A = pd.DataFrame(arows)
    T[T.set.str.endswith("E")].to_csv(f"{OUT}/D_tile_E.csv", index=False); T[~T.set.str.endswith("E")].to_csv(f"{OUT}/D_tile_U.csv", index=False)
    A.to_csv(f"{OUT}/D_A_tile.csv", index=False)
    T = pd.concat([pd.read_csv(f"{OUT}/D_tile_E.csv"), pd.read_csv(f"{OUT}/D_tile_U.csv")], ignore_index=True)
    I = {n: T[(T.set == n) & (T.included == True)].reset_index(drop=True) for n in SETS}
    S = dict(contract_sha256="257ab6ea6d6132429cbd440e063711eb50ff1b4740ac25843189e827e3978af1", smoke=bool(SMOKE),
             n_included={n: int(len(v)) for n, v in I.items()})
    # ---------- D-1
    D1 = {}
    for n in ("O'E", "N'E"):
        D1[n] = dict(M1_E=m1(I[n]), M2_E=m2(I[n]), metrics=block(I[n], A, n, ("T05", "T5", "T10")))
    trig = {n: abs(D1[n]["M2_E"]["S0"]) >= 0.2 for n in D1}
    ok = {n: bool(abs(D1[n]["M2_E"]["S1"]) <= 0.1 and abs(D1[n]["M2_E"]["S1"]) <= abs(D1[n]["M2_E"]["S0"]) / 3) for n in D1}
    if all(trig.values()): h7 = dict(status="testable", holds=bool(all(ok.values())))
    elif not any(trig.values()): h7 = dict(status="untestable", holds=None)
    else: k = [n for n in trig if trig[n]][0]; h7 = dict(status="one_set_only", set=k, holds=ok[k])
    h7.update(triggered=trig, per_set_condition=ok); D1["H7"] = h7
    # U values on the same tiles (Confirm02 / Confirm03 A-C outputs), for side-by-side reporting only
    try:
        U2 = json.load(open(f"{C02}/results/mech/summary.json"))["aggregates"]; U3 = json.load(open(os.path.join(HERE, "..", "results", "summary.json")))
        D1["U_same_tiles"] = dict(aggregates={k: v for k, v in U2.items() if k.split("|")[0] in ("O'", "N'") and k.split("|")[1] in ("C0", "C1", "C2")},
                                  M2_U=json.load(open(f"{C02}/results/mech/summary.json"))["M2"],
                                  H123={k: U3[k] for k in ("H1", "H2", "H3") if k in U3})
    except Exception as e: D1["U_same_tiles"] = dict(error=type(e).__name__)
    S["D1"] = D1
    # ---------- D-2
    U = pd.concat([I["N'"], I["S''"]], ignore_index=True)
    H8 = {f"{mode}_{var}": f2(U, f"{mode}__mstar_{var}") for mode in ("C1", "C0") for var in ("cw", "pub")}
    H8["judged_on"] = "C1_cw"; H8["holds"] = H8["C1_cw"]["pass"]
    H8["ratio_ge_3_evaluated_and_met"] = bool(H8["C1_cw"]["ratio_evaluated"] and (H8["C1_cw"]["ratio"] or 0) >= 3)
    try:
        P1 = pd.read_csv(f"{RUN}/part1_Nprime/by_tile.csv"); P1 = P1[P1.included == True][["tile", "m"]].merge(I["N'"][["tile", "C0__mstar_cw"]], on="tile")
        H8["Nprime_sample_m_vs_fulltile_m_spearman"] = float(spearmanr(P1.m, P1.C0__mstar_cw).statistic); H8["Nprime_n_compared"] = int(len(P1))
    except Exception as e: H8["Nprime_sample_m_vs_fulltile_m_spearman"] = type(e).__name__
    Np_sp = I["N'"][I["N'"].dclass == "sparse"].reset_index(drop=True)
    S["D2"] = dict(H8=H8, Spp=dict(M1=m1(I["S''"]), M2=m2(I["S''"]), n_by_class={k: int(v) for k, v in I["S''"].dclass.value_counts().items()},
                                   metrics=block(I["S''"], A, "S''", ("T05", "T5", "T10"))),
                   Nprime_sparse=dict(n_tiles=int(len(Np_sp)), metrics=block(Np_sp, A, "N'", ("T05", "T5", "T10")) if len(Np_sp) >= 2 else None))
    # ---------- D0 reuse checks on N' (U): same numbers as Confirm03 s96 A_tile and Confirm02 by_tile_N' (booleans only)
    try:
        a3 = pd.read_csv(os.path.join(HERE, "..", "results", "A_tile.csv")); a3 = a3[a3.set == "N'"]
        kc = ["tile", "var", "mode", "T", "dir"]; vc = ["n_mov", "n_lost", "E1", "E2", "n_pool"]
        mg = A[A.set == "N'"][kc + vc].merge(a3[kc + vc], on=kc, suffixes=("", "_ref"))
        S["D0_partA_equal_s96"] = dict(n_rows=int(len(mg)), equal=bool(len(mg) > 0 and all(np.allclose(mg[c], mg[c + "_ref"], rtol=0, atol=1e-9) for c in vc)))
        b2 = pd.read_csv(f"{C02}/results/mech/by_tile_N'.csv"); b2 = b2[b2.included == True].set_index("tile")
        cols = [c for c in I["N'"].columns if re.match(r"C[012]__(pub|cw)_T(5|10|05)_(dis|loss|gain|persA|built_dis|built_flag)$", c) and c in b2.columns]
        J = I["N'"].set_index("tile"); J = J[J.index.isin(b2.index)]
        S["D0_metrics_equal_s95"] = dict(n_tiles=int(len(J)), n_cols=len(cols),
                                       equal=bool(len(J) > 0 and all(np.allclose(J[c].values, b2.loc[J.index, c].values, rtol=0, atol=1e-9, equal_nan=True) for c in cols)))
    except Exception as e: S["D0_error"] = f"{type(e).__name__}: {e}"
    json.dump(S, open(f"{OUT}/summary_D.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print("written", OUT, flush=True)
