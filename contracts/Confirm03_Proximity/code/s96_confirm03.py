"""Confirm03_Proximity s96 (cloud): parts A, B, C of CONTRACT_CN.md (SHA 257ab6ea...8af1). Part D runs later (s97/s98).
Sets and inclusion as Confirm01/02: O (P1, dg = 0), O' and N' (P2, dg = gB - gA from Confirm02 s94).
Variants: C0 (vB as published) for all sets, C1 (vB - dg) for O' and N'; velocities pub (primary) and cw (secondary).
Arithmetic mirrors Confirm02 (float64 from float32 arrays, C1 = vB - dg via rsc_core.harmonize, no rounding), so C0/C1 losses equal
the Confirm01/02 tile values (check A0).
Part A, per tile, threshold T in {5, 10}, direction 'loss' (reference A, other B') and 'gain' (reference B', other A):
  moving = |ref| >= T; lost = |other| < T; d = |ref| - T; a = |other - ref|;
  bins of d: [0,0.5) [0.5,1) [1,2) [2,5) [5,inf);
  N1: E1 = 0.5 * 1(| |ref| - a | < T);
  N2: eps pool = (other - ref) over points of the tile with |ref| < T/2; E2 = #{eps : -T-ref < eps < T-ref} / #pool, exact by sorted search.
  Interpretation fixed here before any output was read: a tile whose pool has < 50 points contributes nothing to the N2 sums
  (its moving points are left out of L_obs used for rho2 as well); counts are reported.
  Aggregates: pooled sums over tiles; ratios of sums; tile bootstrap 1000, default_rng(20261120), 2.5/97.5 percentiles.
Part B (stations, common-window station velocities as Confirm01 F8): B-1 tile medians on all points (pub);
  B-2 per-release offsets c_r = median_s(v_r,s - v_GNSS,s); dg_s = median of point dg within the station's matched points;
  station bootstrap 1000, default_rng(20261121). B-3 decision-level R at stations (T0, T0.5), EGMS minus common offset c.
Part C: N' strata (manifest stratum), report only; T10 built-up s_g from Confirm02 by_tile_N'.csv (C0, C1 columns) with
  tile bootstrap default_rng(20261120).
Usage: python s96_confirm03.py OUT C01_RUN C02_DIR DATA_DIR [NPROC] [SMOKE_TILES comma-separated: tiles and stations restricted to them]"""
import sys, os, glob, json, re
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import theilslopes
OUT, RUN, C02, DATA = sys.argv[1:5]; NPROC = int(sys.argv[5]) if len(sys.argv) > 5 else 2
SMOKE = sys.argv[6].split(",") if len(sys.argv) > 6 else None
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, f"{C02}/code/rsc"); import rsc_core as rc
SETS = {"O": (f"{RUN}/fulltile_O", None), "O'": (f"{RUN}/fulltile_Oprime", f"{C02}/data/g_Oprime"), "N'": (f"{RUN}/fulltile_Nprime", f"{C02}/data/g_Nprime")}
SUB = {"O'": "Oprime", "N'": "Nprime"}
BINS = [0, 0.5, 1, 2, 5, np.inf]; NB = len(BINS) - 1; NBOOT = 1000; MINPOOL = 50
def included(m): return ("excluded" not in m) and m["nA_epochs"] >= 60 and m["nB_epochs"] >= 60 and m["W_days"] >= 1095 and m["n_matched"] >= 50
def load_dg(name, t, n):
    if SETS[name][1] is None: return np.zeros(n)
    z = np.load(f"{SETS[name][1]}/{t}.npz"); return z["gB"].astype(np.float64) - z["gA"].astype(np.float64)

def partA(ref, oth, T):
    f = np.isfinite(ref) & np.isfinite(oth); ref, oth = ref[f], oth[f]
    mov = np.abs(ref) >= T; r = ref[mov]; o = oth[mov]; lost = np.abs(o) < T; d = np.abs(r) - T; a = np.abs(o - r)
    out = dict(n_mov=int(mov.sum()), n_lost=int(lost.sum()))
    b = np.digitize(d, BINS[1:-1])
    for k in range(NB): out[f"mov_b{k}"] = int((b == k).sum()); out[f"lost_b{k}"] = int((lost & (b == k)).sum())
    out["E1"] = float(0.5 * (np.abs(np.abs(r) - a) < T).sum())
    pool = np.sort((oth - ref)[np.abs(ref) < T / 2]); out["n_pool"] = int(len(pool))
    if len(pool) >= MINPOOL and len(r):
        hi = np.searchsorted(pool, T - r, side="left"); lo = np.searchsorted(pool, -T - r, side="right")
        out["E2"] = float((np.clip(hi - lo, 0, None) / len(pool)).sum()); out["n_mov_N2"] = out["n_mov"]; out["n_lost_N2"] = out["n_lost"]
    else:
        out["E2"] = 0.0; out["n_mov_N2"] = 0; out["n_lost_N2"] = 0
    return out

def work(args):
    name, t = args; fd = SETS[name][0]; m = json.load(open(f"{fd}/{t}.json"))
    if not included(m): return [], dict(tile=t, set=name, included=False)
    z = np.load(f"{fd}/{t}.npz"); E, N = z["E"], z["N"]; dg = load_dg(name, t, len(E)); rows = []
    vp = z["vB_pub"].astype(float) - z["vA_pub"].astype(float); fp = np.isfinite(vp)
    b1 = dict(tile=t, set=name, included=True, n=min(m["nA_epochs"], m["nB_epochs"]), med_dv=float(np.median(vp[fp])),
              med_dg=float(np.median(dg)), med_dv_minus_dg=float(np.median((vp - dg)[fp])))
    for var in ("pub", "cw"):
        vA = z[f"vA_{var}"].astype(float); vB = z[f"vB_{var}"].astype(float)
        for mode in (["C0"] if name == "O" else ["C0", "C1"]):
            vBh = rc.harmonize(E, N, vA, vB, dg, mode)
            for T in (5.0, 10.0):
                for dirn, ref, oth in (("loss", vA, vBh), ("gain", vBh, vA)):
                    r = partA(ref, oth, T); r.update(set=name, tile=t, var=var, mode=mode, T=int(T), dir=dirn); rows.append(r)
    return rows, b1

def ratio_ci(df, num, den, seed):
    tiles = df.tile.unique(); g = np.random.default_rng(seed); S = df.groupby("tile")[[num, den]].sum().reindex(tiles)
    a, b = S[num].values, S[den].values; est = a.sum() / b.sum() if b.sum() > 0 else np.nan; bs = []
    for _ in range(NBOOT):
        ix = g.integers(0, len(tiles), len(tiles)); bs.append(a[ix].sum() / b[ix].sum() if b[ix].sum() > 0 else np.nan)
    return float(est), [float(np.nanpercentile(bs, 2.5)), float(np.nanpercentile(bs, 97.5))]

def summarize(df):
    """df: rows of one set/var/mode/T/dir."""
    s = df.sum(numeric_only=True); o = dict(n_tiles=int(df.tile.nunique()), n_mov=int(s.n_mov), n_lost=int(s.n_lost))
    o["L_obs"] = s.n_lost / s.n_mov if s.n_mov else np.nan; o["L1"] = s.E1 / s.n_mov if s.n_mov else np.nan
    o["L2"] = s.E2 / s.n_mov_N2 if s.n_mov_N2 else np.nan; o["n_mov_N2"] = int(s.n_mov_N2)
    d = df.assign(E1r=df.E1); o["rho1"], o["rho1_ci"] = ratio_ci(d, "n_lost", "E1", 20261120)
    d2 = df[df.n_mov_N2 > 0]
    if len(d2) and d2.E2.sum() > 0: o["rho2"], o["rho2_ci"] = ratio_ci(d2, "n_lost_N2", "E2", 20261120)
    else: o["rho2"], o["rho2_ci"] = np.nan, [np.nan, np.nan]
    o["bins"] = []
    for k in range(NB):
        mv, ls = int(s[f"mov_b{k}"]), int(s[f"lost_b{k}"])
        o["bins"].append(dict(bin=f"[{BINS[k]},{BINS[k+1]})", n_mov=mv, n_lost=ls, rate=ls / mv if mv else np.nan, share_of_losses=ls / s.n_lost if s.n_lost else np.nan))
    o["share_lost_d_lt_1"] = (s.lost_b0 + s.lost_b1) / s.n_lost if s.n_lost else np.nan
    o["rate_b0_over_rate_b3"] = (o["bins"][0]["rate"] / o["bins"][3]["rate"]) if (o["bins"][3]["rate"] and o["bins"][3]["rate"] > 0) else np.nan
    return o

if __name__ == "__main__":
    A_rows, B1 = [], []
    for name, (fd, gd) in SETS.items():
        tiles = sorted(os.path.basename(p)[:-5] for p in glob.glob(f"{fd}/*.json"))
        if SMOKE: tiles = [t for t in tiles if t in SMOKE]
        with Pool(NPROC) as p: res = p.map(work, [(name, t) for t in tiles], chunksize=1)
        for rows, b in res: A_rows += rows; B1.append(b)
        print(name, "tiles", len(tiles), flush=True)
    A = pd.DataFrame(A_rows); A.to_csv(f"{OUT}/A_tile.csv", index=False); B1 = pd.DataFrame(B1); B1.to_csv(f"{OUT}/B1_tile.csv", index=False)
    S = dict(contract_sha256="257ab6ea6d6132429cbd440e063711eb50ff1b4740ac25843189e827e3978af1")
    # ---- A0 check: pooled loss/gain counts equal to the Confirm02 tile tables (pub, C0/C1)
    chk = {}
    for name in SETS:
        bt = pd.read_csv(f"{C02}/results/mech/by_tile_{name}.csv"); bt = bt[bt.included == True]
        for mode in (["C0"] if name == "O" else ["C0", "C1"]):
            for T in (5, 10):
                sub = A[(A.set == name) & (A["mode"] == mode) & (A["var"] == "pub") & (A["T"] == T) & (A.dir == "loss")].set_index("tile")
                col_l, col_n = f"{mode}__pub_T{T}_loss", f"{mode}__pub_T{T}_n_movA"
                if col_l not in bt.columns: continue
                b = bt.set_index("tile").reindex(sub.index); lostB = (b[col_l].fillna(0) * b[col_n]).round().astype(int)
                chk[f"{name}|{mode}|T{T}"] = bool((lostB.values == sub.n_lost.values).all() and (b[col_n].values == sub.n_mov.values).all())
    S["A0_equal_to_Confirm02"] = chk
    # ---- A aggregates
    AG = {}; binrows = []
    for (name, var, mode, T, dirn), g in A.groupby(["set", "var", "mode", "T", "dir"]):
        o = summarize(g); AG[f"{name}|{var}|{mode}|T{T}|{dirn}"] = o
        for b in o["bins"]: binrows.append(dict(set=name, var=var, mode=mode, T=T, dir=dirn, **b))
    pd.DataFrame(binrows).to_csv(f"{OUT}/A_bins.csv", index=False); S["A"] = AG
    J = {"O": "C0", "O'": "C1", "N'": "C1"}; key = lambda n, d: f"{n}|pub|{J[n]}|T5|{d}"
    h1 = {n: dict(share=AG[key(n, "loss")]["share_lost_d_lt_1"], ratio=AG[key(n, "loss")]["rate_b0_over_rate_b3"]) for n in J}
    S["H1"] = dict(per_set=h1, holds=bool(all(v["share"] >= 0.5 and v["ratio"] >= 3 for v in h1.values())))
    def lab(o):
        lo, hi = o["rho1_ci"]
        if lo > 1: return "above"
        if hi < 1: return "below"
        return "includes_1"
    h2 = {}
    for n in J:
        L, G = AG[key(n, "loss")], AG[key(n, "gain")]
        both_in = all(0.8 <= o["rho1"] <= 1.25 and o["rho1_ci"][0] <= 1 <= o["rho1_ci"][1] for o in (L, G))
        verdict = "both_above" if (lab(L) == "above" and lab(G) == "above") else ("consistent_with_random_sign" if both_in else ("directional" if (lab(L) == "above") != (lab(G) == "above") else "other"))
        h2[n] = dict(rho1_loss=L["rho1"], ci_loss=L["rho1_ci"], rho1_gain=G["rho1"], ci_gain=G["rho1_ci"], verdict=verdict)
    S["H2"] = dict(per_set=h2, holds=bool(all(v["verdict"] == "both_above" for v in h2.values())))
    h3 = {n: dict(rho2=AG[key(n, "loss")]["rho2"], ci=AG[key(n, "loss")]["rho2_ci"]) for n in J}
    S["H3"] = dict(per_set=h3, holds=bool(all(v["rho2"] >= 1.5 and v["ci"][0] > 1 for v in h3.values())))
    # ---- B-1
    b1 = {}
    for name in ("O'", "N'"):
        d = B1[(B1.set == name) & (B1.included == True)]
        on = lambda x: float(np.mean(np.abs(x * 10 - np.round(x * 10)) < 1e-6))
        sl = theilslopes(d.med_dv, d.med_dg)
        b1[name] = dict(n_tiles=int(len(d)), median_of_med_dv=float(d.med_dv.median()), median_of_med_dg=float(d.med_dg.median()),
                        median_of_med_dv_minus_dg=float(d.med_dv_minus_dg.median()), frac_abs_med_dv_minus_dg_le_0p1=float((d.med_dv_minus_dg.abs() <= 0.1 + 1e-9).mean()),
                        theil_sen_slope=float(sl[0]), theil_sen_ci=[float(sl[2]), float(sl[3])],
                        frac_on_0p1_grid=dict(med_dv=on(d.med_dv.values), med_dg=on(d.med_dg.values), med_dv_minus_dg=on(d.med_dv_minus_dg.values)))
    S["B1"] = b1
    # ---- B-2 / B-3 stations
    sta = {}
    for name in ("O'", "N'"):
        P = pd.read_csv(f"{RUN}/part1_{SUB[name]}/by_point.csv"); P = P[P.kind == "sta"].copy(); P["dg"] = np.nan
        if SMOKE: P = P[P.tile.isin(SMOKE)]
        for t, idx in P.groupby("tile").groups.items():
            z = np.load(f"{SETS[name][0]}/{t}.npz"); dg = load_dg(name, t, len(z["E"])); kf = rc.key(z["E"], z["N"]); o = np.argsort(kf)
            kp = rc.key(P.loc[idx, "E"].values, P.loc[idx, "N"].values); i = np.clip(np.searchsorted(kf[o], kp), 0, len(kf) - 1); ok = kf[o][i] == kp
            P.loc[idx, "dg"] = np.where(ok, dg[o][i], np.nan)
        dgs = P.groupby("sta").dg.median().rename("dg_s"); unm = int(P.dg.isna().sum())
        ST = pd.read_csv(f"{RUN}/part1_{SUB[name]}/by_station.csv")[["sta", "tile", "dclass", "vA", "vB", "v_gnss"]].merge(dgs, on="sta", how="left")
        if SMOKE: ST = ST[ST.tile.isin(SMOKE)]
        sta[name] = (ST, unm)
    O_tiles = set(pd.read_csv(f"{DATA}/manifest_O.csv").tile)
    P1 = pd.read_csv(f"{DATA}/relcons_by_station.csv"); P1 = P1[P1.tile.isin(O_tiles)][["sta", "tile", "dclass", "vA", "vB", "v_gnss"]]
    if SMOKE: P1 = P1[P1.tile.isin(SMOKE)]
    def boot_st(n, seed=20261121):
        g = np.random.default_rng(seed); return [g.integers(0, n, n) for _ in range(NBOOT)]
    def ci(vals): return [float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))]
    def cr(v, vg): return float(np.median(v - vg))
    B2 = {}; strows = []
    for name in ("O'", "N'"):
        ST, unm = sta[name]; ST = ST[ST.dg_s.notna()].reset_index(drop=True); bi = boot_st(len(ST))
        def f(D):
            cA = cr(D.vA, D.v_gnss); cB = cr(D.vB, D.v_gnss); cB1 = cr(D.vB - D.dg_s, D.v_gnss); mdg = float(np.median(D.dg_s))
            return cA, cB, cB1, mdg, (cB - cA) - mdg, cB1 - cA
        est = f(ST); bs = np.array([f(ST.iloc[ix]) for ix in bi])
        names = ["c_2020_24_P2", "c_2020_25", "c_2020_25_C1", "median_dg_s", "x_H4a", "x_H4b"]
        B2[name] = {k: dict(est=float(e), ci=ci(bs[:, j])) for j, (k, e) in enumerate(zip(names, est))}
        B2[name].update(n_stations=int(len(ST)), n_points_unmatched=unm,
                        H4=bool(abs(est[4]) <= 0.20 and ci(bs[:, 4])[0] <= 0 <= ci(bs[:, 4])[1] and abs(est[5]) <= 0.20))
        strows.append(ST.assign(set=name))
    bi = boot_st(len(P1))
    f1 = lambda D: (cr(D.vA, D.v_gnss), cr(D.vB, D.v_gnss), cr(D.vB, D.v_gnss) - cr(D.vA, D.v_gnss))
    est = f1(P1); bs = np.array([f1(P1.iloc[ix]) for ix in bi])
    B2["P1"] = dict(n_stations=int(len(P1)), c_2019_23=dict(est=est[0], ci=ci(bs[:, 0])), c_2020_24_P1=dict(est=est[1], ci=ci(bs[:, 1])),
                    diff=dict(est=est[2], ci=ci(bs[:, 2])), c_mid=float(np.median((P1.vA + P1.vB) / 2 - P1.v_gnss)), H5=bool(abs(est[2]) <= 0.20))
    STo = sta["O'"][0]; common = sorted(set(P1.sta) & set(STo.sta)); a = P1.set_index("sta").loc[common]; b = STo.set_index("sta").loc[common]
    x = cr(a.vB, a.v_gnss) - cr(b.vA, b.v_gnss); bi = boot_st(len(common))
    bs = [cr(a.vB.iloc[ix], a.v_gnss.iloc[ix]) - cr(b.vA.iloc[ix], b.v_gnss.iloc[ix]) for ix in bi]
    dgc = b.dg_s.dropna()
    B2["H6"] = dict(n_common=len(common), diff_2020_24_P1_minus_P2=float(x), ci=ci(bs), holds=bool(abs(x) <= 0.20),
                    report_P1_mid_offset_common=float(np.median((a.vA + a.vB) / 2 - a.v_gnss)), report_median_dg_s_common=float(dgc.median()) if len(dgc) else None,
                    report_implied_ETRF2020_offset_of_2019_23=float(cr(a.vA, a.v_gnss) + dgc.median()) if len(dgc) else None)
    S["B2"] = B2
    pd.concat(strows + [P1.assign(set="O_P1")]).to_csv(f"{OUT}/B_stations.csv", index=False)
    def cls(v, th): return (v < 0).astype(int) if th == "T0" else np.where(v < -0.5, -1, np.where(v > 0.5, 1, 0))
    def rdec(D, vBcol, th):
        c = float(np.median((D.vA + D[vBcol]) / 2 - D.v_gnss)); a_ = cls(D.vA - c, th); b_ = cls(D[vBcol] - c, th); g_ = cls(D.v_gnss, th)
        den = 0.5 * (np.mean(a_ != g_) + np.mean(b_ != g_)); return float(np.mean(a_ != b_) / den) if den > 0 else np.nan, float(np.mean(a_ != b_)), float(den)
    B3 = {}
    for name in ("O'", "N'"):
        ST = sta[name][0]; ST = ST[ST.dg_s.notna()].assign(vB1=lambda D: D.vB - D.dg_s).reset_index(drop=True); bi = boot_st(len(ST))
        for th in ("T0", "T05"):
            for vc, md in (("vB", "C0"), ("vB1", "C1")):
                e = rdec(ST, vc, th); bs = [rdec(ST.iloc[ix], vc, th)[0] for ix in bi]
                B3[f"{name}|{md}|{th}"] = dict(R_dec=e[0], ci=ci(bs), p_release_flip=e[1], p_egms_gnss_flip=e[2])
    bi = boot_st(len(P1))
    for th in ("T0", "T05"):
        e = rdec(P1, "vB", th); bs = [rdec(P1.iloc[ix], "vB", th)[0] for ix in bi]
        B3[f"O_P1|C0|{th}"] = dict(R_dec=e[0], ci=ci(bs), p_release_flip=e[1], p_egms_gnss_flip=e[2])
    S["B3"] = B3
    # ---- C N' strata
    man = pd.read_csv(f"{DATA}/manifest_Nprime.csv")[["tile", "stratum"]]; bt = pd.read_csv(f"{C02}/results/mech/by_tile_N'.csv"); bt = bt[bt.included == True].merge(man, on="tile")
    AN = A[A.set == "N'"].merge(man, on="tile"); C = {}; crow = []
    for st_, g in AN.groupby("stratum"):
        o = dict(n_tiles=int(g.tile.nunique()))
        for mode in ("C0", "C1"):
            for T in (5, 10):
                sub = g[(g["var"] == "pub") & (g["mode"] == mode) & (g["T"] == T) & (g.dir == "loss")]
                o[f"loss_pooled_{mode}_T{T}"] = float(sub.n_lost.sum() / sub.n_mov.sum()) if sub.n_mov.sum() else np.nan
        b = bt[bt.stratum == st_]
        for mode in ("C0", "C1"):
            o[f"A3_T10_{mode}"] = float(b[f"{mode}__pub_T10_built_dis"].sum() / b[f"{mode}__pub_T10_built_flag"].sum()) if b[f"{mode}__pub_T10_built_flag"].sum() > 0 else np.nan
        bb = b[(b["C0__pub_T10_built_flag"] > 0) & (b["C0__pub_T10_built_dis"] > 0)].reset_index(drop=True)
        def sg(J):
            d0 = J["C0__pub_T10_built_dis"].sum() / J["C0__pub_T10_built_flag"].sum(); d1 = J["C1__pub_T10_built_dis"].sum() / J["C1__pub_T10_built_flag"].sum(); return (d0 - d1) / d0
        if len(bb):
            gg = np.random.default_rng(20261120); bs = [sg(bb.iloc[gg.integers(0, len(bb), len(bb))]) for _ in range(NBOOT)]
            o["sg_A3_T10"] = float(sg(bb)); o["sg_A3_T10_ci"] = ci(bs); o["sg_n_tiles"] = int(len(bb))
        sub = g[(g["var"] == "pub") & (g["mode"] == "C1") & (g["T"] == 5)]
        if sub[sub.dir == "loss"].n_lost.sum() > 0:
            Ls = summarize(sub[sub.dir == "loss"]); Gs = summarize(sub[sub.dir == "gain"])
            o.update(H1_share_d_lt_1=Ls["share_lost_d_lt_1"], H1_rate_ratio=Ls["rate_b0_over_rate_b3"], H2_rho1_loss=Ls["rho1"], H2_rho1_loss_ci=Ls["rho1_ci"],
                     H2_rho1_gain=Gs["rho1"], H2_rho1_gain_ci=Gs["rho1_ci"], H3_rho2=Ls["rho2"], H3_rho2_ci=Ls["rho2_ci"])
        C[st_] = o; crow.append(dict(stratum=st_, **{k: v for k, v in o.items() if not isinstance(v, list)}))
    S["C"] = C; pd.DataFrame(crow).to_csv(f"{OUT}/C_strata.csv", index=False)
    json.dump(S, open(f"{OUT}/summary.json", "w"), indent=1, default=float); print("done", flush=True)
