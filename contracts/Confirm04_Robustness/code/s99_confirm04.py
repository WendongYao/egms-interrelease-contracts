"""Confirm04_Robustness s99 (cloud): CONTRACT_CN.md (SHA c57fbc7e...0a58). R1 rounding sensitivity, R2 spatial-block bootstrap,
R3 threshold curves, R4 uncertainty buffer (exploratory). No earlier verdict is changed here; outputs only feed the wording rules.
Reused definitions: Confirm03 s96 partA (imported); Confirm02 s95 aggD (exec of its verbatim source); rsc_core.harmonize; tile
inclusion as Confirm02/03 (not excluded, both releases >= 60 epochs in W, W >= 1095 d, >= 50 matched points).
Interpretations fixed here before any output was read:
  R1  Harrell-Davis median: weights w_i = I_{i/n}(a,a) - I_{(i-1)/n}(a,a), a = (n+1)/2, on sorted finite values; used both within the
      tile and across tiles. Dither: for every point, v_A, v_B, g_A, g_B get independent U(-0.05, 0.05); generator per tile
      default_rng([20261123, set index, tile index in the sorted tile list]); 20 replicates, replicate r uses draws r of each tile.
      Tile medians use points with finite published dv (as Confirm02 s95). Tile-mean variant: mean within tile, median across tiles.
  R2  Blocks (E//3, N//3) from the tile name (600 km: E//6, N//6); stations take the block of their tile. One generator
      default_rng(20261122) per statistic; each replicate draws, with replacement, as many blocks as the statistic's data contain and
      takes every tile (station) of the drawn blocks. Statistics are the original definitions: rho1/rho2 = ratio of pooled sums
      (A_tile.csv); x_H4a = (c_B - c_A) - median dg_s on stations; Spearman rho(n, m) for H8 (D_tile_U, C1 cw) and F2 (Confirm01
      part1_Nprime/by_tile, m); s_g via s95 aggD on the s95 share subset (tiles with finite, positive C0 tile value); R_dec as s96.
      A verdict is re-evaluated by substituting the block interval into the original rule, with the original point estimates.
  R3  Tile median loss over tiles with at least one moving point. rho2 pooled over tiles whose noise pool has >= 50 points (s96).
  R4  Moving points of release A whose mean_velocity_std_A is not finite are counted separately and are not put in the buffer.
Outputs: OUT/summary.json, R1.csv, R2.csv, R3_curves.csv, R4.csv. Prints counts only.
Usage: python s99_confirm04.py OUT C01_RUN C02_DIR C03_DIR [NPROC] [SMOKE]"""
import sys, os, glob, json, re, ast
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.special import betainc
from scipy.stats import spearmanr
OUT, RUN, C02, C03 = sys.argv[1:5]; NPROC = int(sys.argv[5]) if len(sys.argv) > 5 else 2; SMOKE = len(sys.argv) > 6
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, f"{C02}/code/rsc"); sys.path.insert(0, f"{C03}/code")
import rsc_core as rc
_a = sys.argv; sys.argv = ["s96", os.path.join(OUT, "_s96_unused"), RUN, C02, f"{C03}/data"]; import s96_confirm03 as s96; sys.argv = _a
SEED, NBOOT = 20261118, 1000
for node in ast.parse(open(f"{C02}/code/s95_mech.py").read()).body:
    if isinstance(node, ast.FunctionDef) and node.name in ("tileD", "aggD"):
        exec(compile(ast.Module([node], []), "s95_mech.py", "exec"), globals())
DD = f"{C03}/dataD"
SETS = {"O'U": (f"{RUN}/fulltile_Oprime", f"{C02}/data/g_Oprime"), "N'U": (f"{RUN}/fulltile_Nprime", f"{C02}/data/g_Nprime"),
        "O'E": (f"{DD}/fulltile_E_Oprime", f"{DD}/gE_Oprime"), "N'E": (f"{DD}/fulltile_E_Nprime", f"{DD}/gE_Nprime"),
        "S''": (f"{DD}/fulltile_Spp", f"{DD}/g_Spp"), "O": (f"{RUN}/fulltile_O", None)}
R1SETS = ["O'U", "N'U", "O'E", "N'E", "S''"]
R3SETS = {"O": "C0", "O'U": "C1", "N'U": "C1", "O'E": "C1", "N'E": "C1"}
TS = [float(t) for t in range(2, 16)]; NDITH = 20
def included(m): return ("excluded" not in m) and m["nA_epochs"] >= 60 and m["nB_epochs"] >= 60 and m["W_days"] >= 1095 and m["n_matched"] >= 50
def tiles_of(name):
    fd = SETS[name][0]; t = sorted(os.path.basename(p)[:-5] for p in glob.glob(f"{fd}/*.json"))
    return t[:2] if SMOKE else t
def hd_median(x):
    x = np.sort(np.asarray(x, float)[np.isfinite(x)]); n = len(x)
    if n == 0: return np.nan
    a = (n + 1) / 2; u = np.arange(n + 1) / n; w = np.diff(betainc(a, a, u)); return float(w @ x)

def work(args):
    name, si, ti, t = args; fd, gd = SETS[name]; m = json.load(open(f"{fd}/{t}.json"))
    if not included(m): return dict(set=name, tile=t, included=False), [], []
    z = np.load(f"{fd}/{t}.npz"); E, N = z["E"], z["N"]
    vA = z["vA_pub"].astype(np.float64); vB = z["vB_pub"].astype(np.float64)
    if gd is None: gA = gB = np.zeros(len(E))
    else: g = np.load(f"{gd}/{t}.npz"); gA = g["gA"].astype(np.float64); gB = g["gB"].astype(np.float64)
    dg = gB - gA; out = dict(set=name, tile=t, included=True)
    if name in R1SETS:
        dv = vB - vA; f = np.isfinite(dv); x0, x1 = dv[f], (dv - dg)[f]
        out.update(med0=float(np.median(x0)), med1=float(np.median(x1)), hd0=hd_median(x0), hd1=hd_median(x1), mean0=float(np.mean(x0)), mean1=float(np.nanmean(x1)))
        rng = np.random.default_rng([20261123, si, ti]); n = len(E)
        for r in range(NDITH):
            uA, uB, ugA, ugB = (rng.uniform(-0.05, 0.05, n) for _ in range(4))
            d = (vB + uB) - (vA + uA); dgr = (gB + ugB) - (gA + ugA); ff = np.isfinite(d)
            out[f"dith0_{r}"] = float(np.median(d[ff])); out[f"dith1_{r}"] = float(np.median((d - dgr)[ff]))
    r3, r4 = [], []
    if name in R3SETS:
        vBh = rc.harmonize(E, N, vA, vB, dg, R3SETS[name]); mvs = z["mvsA"].astype(np.float64)
        for T in TS:
            for dirn, ref, oth in (("loss", vA, vBh), ("gain", vBh, vA)):
                r = s96.partA(ref, oth, T); r.update(set=name, tile=t, T=T, dir=dirn, mode=R3SETS[name])
                r3.append({k: v for k, v in r.items() if not k.startswith(("mov_b", "lost_b"))})
            if T in (5.0, 10.0) and name in ("O", "O'U", "N'U"):
                ok = np.isfinite(vA) & np.isfinite(vBh); a, b, s = vA[ok], vBh[ok], mvs[ok]
                mov = np.abs(a) >= T; lost = mov & (np.abs(b) < T); d = np.abs(a) - T
                fin = np.isfinite(s) & (s > 0) if True else None
                buf = mov & fin & (d < 2 * np.where(fin, s, 0)); b1 = mov & (d < 1)
                r4.append(dict(set=name, tile=t, T=T, n_mov=int(mov.sum()), n_lost=int(lost.sum()), n_mov_mvs_missing=int((mov & ~fin).sum()),
                               n_buf=int(buf.sum()), n_lost_buf=int((lost & buf).sum()), n_d1=int(b1.sum()), n_lost_d1=int((lost & b1).sum())))
    return out, r3, r4

# ---------------- R2 helpers
def blk(tile, k): e, n = [int(x) for x in re.findall(r"\d+", tile)]; return (e // k, n // k)
def block_draws(blocks, seed=20261122, nboot=NBOOT):
    """blocks: array of block labels per row. Returns list of row-index arrays."""
    lab = pd.Series(list(map(str, blocks))); groups = [np.where(lab.values == u)[0] for u in sorted(lab.unique())]
    g = np.random.default_rng(seed); out = []
    for _ in range(nboot):
        pick = g.integers(0, len(groups), len(groups)); out.append(np.concatenate([groups[i] for i in pick]))
    return out, len(groups)
def pct(v): v = np.asarray(v, float); return [float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))]
def ratio_block(df, num, den, k):
    S = df.groupby("tile")[[num, den]].sum(); idx, nb = block_draws([blk(t, k) for t in S.index]); a, b = S[num].values, S[den].values
    est = a.sum() / b.sum(); bs = [a[i].sum() / b[i].sum() if b[i].sum() > 0 else np.nan for i in idx]; return float(est), pct(bs), nb

if __name__ == "__main__":
    jobs = []
    for si, name in enumerate(SETS):
        for ti, t in enumerate(tiles_of(name)): jobs.append((name, si, ti, t))
    with Pool(NPROC) as p: res = p.map(work, jobs, chunksize=1)
    T1 = pd.DataFrame([r for r, _, _ in res]); R3 = pd.DataFrame([x for _, a, _ in res for x in a]); R4 = pd.DataFrame([x for _, _, b in res for x in b])
    print("tiles", len(T1), "included", int(T1.included.sum()), "R3 rows", len(R3), "R4 rows", len(R4), flush=True)
    S = dict(contract_sha256="c57fbc7ed6e4aa01f11cbbc65ddbb5cde83af61fc573a712b01a876ebb640a58", smoke=SMOKE)
    # ---------------- R1
    r1 = {}; rows1 = []
    m2 = lambda s0, s1: bool(abs(s1) <= 0.1 and abs(s1) <= abs(s0) / 3)
    for name in R1SETS:
        J = T1[(T1.set == name) & (T1.included == True)]
        o = dict(n_tiles=int(len(J)), orig=dict(S0=float(np.median(J.med0)), S1=float(np.median(J.med1))),
                 hd=dict(S0=hd_median(J.hd0.values), S1=hd_median(J.hd1.values)), mean=dict(S0=float(np.median(J.mean0)), S1=float(np.median(J.mean1))))
        d0 = [float(np.median(J[f"dith0_{r}"])) for r in range(NDITH)]; d1 = [float(np.median(J[f"dith1_{r}"])) for r in range(NDITH)]
        o["dither"] = dict(S0=d0, S1=d1, n_meet=int(sum(m2(a, b) for a, b in zip(d0, d1))), S1_min=min(d1), S1_max=max(d1))
        o["hd_meets"] = m2(o["hd"]["S0"], o["hd"]["S1"]); o["dither_meets"] = bool(o["dither"]["n_meet"] >= 19)
        r1[name] = o; rows1.append(dict(set=name, n_tiles=o["n_tiles"], S0=o["orig"]["S0"], S1=o["orig"]["S1"], S0_hd=o["hd"]["S0"], S1_hd=o["hd"]["S1"],
                                        S0_mean=o["mean"]["S0"], S1_mean=o["mean"]["S1"], S1_dith_min=o["dither"]["S1_min"], S1_dith_max=o["dither"]["S1_max"], n_dith_meet=o["dither"]["n_meet"]))
    r1["holds"] = bool(all(r1[n]["hd_meets"] and r1[n]["dither_meets"] for n in R1SETS)); S["R1"] = r1; pd.DataFrame(rows1).to_csv(f"{OUT}/R1.csv", index=False)
    # ---------------- R2
    r2 = {}; rows2 = []
    def rec(key, est, ci3, ci6, nb3, nb6, orig_v, new_v):
        r2[key] = dict(est=est, ci_block300=ci3, ci_block600=ci6, n_blocks300=nb3, n_blocks600=nb6, verdict_orig=orig_v, verdict_block300=new_v, unchanged=(orig_v == new_v) if orig_v is not None else None)
        rows2.append(dict(key=key, est=est, lo300=ci3[0], hi300=ci3[1], lo600=ci6[0], hi600=ci6[1], n_blocks300=nb3, verdict_orig=orig_v, verdict_block300=new_v))
    A = pd.read_csv(f"{C03}/results/A_tile.csv"); S3 = json.load(open(f"{C03}/results/summary.json"))
    J = {"O": "C0", "O'": "C1", "N'": "C1"}
    for n, mode in J.items():
        sub = lambda d: A[(A.set == n) & (A["var"] == "pub") & (A["mode"] == mode) & (A["T"] == 5) & (A.dir == d)]
        for d in ("loss", "gain"):
            e, c3, nb3 = ratio_block(sub(d), "n_lost", "E1", 3); _, c6, nb6 = ratio_block(sub(d), "n_lost", "E1", 6)
            o = S3["A"][f"{n}|pub|{mode}|T5|{d}"]
            rec(f"H2|{n}|{d}", e, c3, c6, nb3, nb6, bool(o["rho1_ci"][0] > 1), bool(c3[0] > 1))
        s2 = sub("loss"); s2 = s2[s2.n_mov_N2 > 0]
        e, c3, nb3 = ratio_block(s2, "n_lost_N2", "E2", 3); _, c6, nb6 = ratio_block(s2, "n_lost_N2", "E2", 6); o = S3["A"][f"{n}|pub|{mode}|T5|loss"]
        rec(f"H3|{n}", e, c3, c6, nb3, nb6, bool(o["rho2"] >= 1.5 and o["rho2_ci"][0] > 1), bool(o["rho2"] >= 1.5 and c3[0] > 1))
    B = pd.read_csv(f"{C03}/results/B_stations.csv")
    cr = lambda v, vg: float(np.median(v - vg))
    for n in ("O'", "N'"):
        ST = B[(B.set == n) & B.dg_s.notna()].reset_index(drop=True)
        xa = lambda D: (cr(D.vB, D.v_gnss) - cr(D.vA, D.v_gnss)) - float(np.median(D.dg_s))
        xb = lambda D: cr(D.vB - D.dg_s, D.v_gnss) - cr(D.vA, D.v_gnss)
        est = xa(ST); cis = {}
        for k in (3, 6):
            idx, nb = block_draws([blk(t, k) for t in ST.tile]); cis[k] = (pct([xa(ST.iloc[i]) for i in idx]), nb)
        rule = lambda ci: bool(abs(est) <= 0.20 and ci[0] <= 0 <= ci[1] and abs(xb(ST)) <= 0.20)
        rec(f"H4|{n}", est, cis[3][0], cis[6][0], cis[3][1], cis[6][1], bool(S3["B2"][n]["H4"]), rule(cis[3][0]))
    def rho_block(df, xc, yc, k):
        idx, nb = block_draws([blk(t, k) for t in df.tile]); bs = [spearmanr(df[xc].iloc[i], df[yc].iloc[i]).statistic if df[xc].iloc[i].nunique() > 1 else np.nan for i in idx]
        return float(spearmanr(df[xc], df[yc]).statistic), pct(bs), nb
    DU = pd.read_csv(f"{C03}/resultsD/D_tile_U.csv"); DU = DU[DU.included == True].reset_index(drop=True); SD = json.load(open(f"{C03}/resultsD/summary_D.json"))["D2"]["H8"]["C1_cw"]
    e, c3, nb3 = rho_block(DU, "n", "C1__mstar_cw", 3); _, c6, nb6 = rho_block(DU, "n", "C1__mstar_cw", 6)
    ratio_ok = (not SD["ratio_evaluated"]) or (SD["ratio"] >= 3)
    rec("H8|N'+S''|C1_cw", e, c3, c6, nb3, nb6, bool(SD["pass"]), bool(e <= -0.5 and c3[1] < -0.2 and ratio_ok))
    P1 = pd.read_csv(f"{RUN}/part1_Nprime/by_tile.csv"); P1 = P1[P1.included == True].reset_index(drop=True); F2o = json.load(open(f"{RUN}/part1_Nprime/summary.json"))["F2"]
    e, c3, nb3 = rho_block(P1, "n", "m", 3); _, c6, nb6 = rho_block(P1, "n", "m", 6)
    rec("F2|N'", e, c3, c6, nb3, nb6, bool(F2o["pass"]), bool(e <= -0.5 and c3[1] < -0.2))
    # M5 a / b (Confirm02), shares via s95 aggD on the s95 share subset
    S2 = json.load(open(f"{C02}/results/mech/summary.json"))
    def sg_block(I, met, th, k):
        d0 = tileD(I, "C0", "pub", met, th); Sb = I[np.isfinite(d0) & (d0 > 0)].reset_index(drop=True)
        f = lambda ix: (lambda D0, D1: (D0 - D1) / D0)(aggD(Sb, "C0", "pub", met, th, ix), aggD(Sb, "C1", "pub", met, th, ix))
        idx, nb = block_draws([blk(t, k) for t in Sb.tile]); return float(f(None)), pct([f(i) for i in idx]), nb
    TB = {n: pd.read_csv(f"{C02}/results/mech/by_tile_{n}.csv") for n in ("O'", "N'")}; TB = {n: v[v.included == True].reset_index(drop=True) for n, v in TB.items()}
    a_new, b_new = True, True
    for n in ("O'", "N'"):
        for th in ("T0", "T05"):
            e, c3, nb3 = sg_block(TB[n], "A1", th, 3); _, c6, nb6 = sg_block(TB[n], "A1", th, 6); o = S2["M5"]["shares"][f"{n}|pub|A1|{th}"]
            v0 = bool(o["s_g"] >= 0.3 and o["s_g_ci"][0] > 0); v1 = bool(o["s_g"] >= 0.3 and c3[0] > 0); a_new &= v1
            rec(f"M5a|{n}|A1|{th}", e, c3, c6, nb3, nb6, v0, v1)
        for met in ("A1", "A2", "A3"):
            e, c3, nb3 = sg_block(TB[n], met, "T10", 3); _, c6, nb6 = sg_block(TB[n], met, "T10", 6); o = S2["M5"]["shares"][f"{n}|pub|{met}|T10"]
            v0 = bool(o["s_g_ci"][1] <= 0.10); v1 = bool(c3[1] <= 0.10); b_new &= v1
            rec(f"M5b|{n}|{met}|T10", e, c3, c6, nb3, nb6, v0, v1)
    r2["M5a_set_verdict"] = dict(orig=bool(S2["M5"]["a_holds"]), block300=bool(a_new)); r2["M5b_set_verdict"] = dict(orig=bool(S2["M5"]["b_holds"]), block300=bool(b_new))
    # report only: R_dec, M2 rho, E H2/H3
    def cls(v, th): return (v < 0).astype(int) if th == "T0" else np.where(v < -0.5, -1, np.where(v > 0.5, 1, 0))
    def rdec(D, vc, th):
        c = float(np.median((D.vA + D[vc]) / 2 - D.v_gnss)); a_ = cls(D.vA - c, th); b_ = cls(D[vc] - c, th); g_ = cls(D.v_gnss, th)
        den = 0.5 * (np.mean(a_ != g_) + np.mean(b_ != g_)); return float(np.mean(a_ != b_) / den) if den > 0 else np.nan
    for n in ("O'", "N'"):
        ST = B[(B.set == n) & B.dg_s.notna()].assign(vB1=lambda D: D.vB - D.dg_s).reset_index(drop=True)
        for th in ("T0", "T05"):
            for vc, md in (("vB", "C0"), ("vB1", "C1")):
                cis = {}
                for k in (3, 6):
                    idx, nb = block_draws([blk(t, k) for t in ST.tile]); cis[k] = (pct([rdec(ST.iloc[i], vc, th) for i in idx]), nb)
                rec(f"Rdec|{n}|{md}|{th}", rdec(ST, vc, th), cis[3][0], cis[6][0], cis[3][1], cis[6][1], None, None)
        e, c3, nb3 = rho_block(TB[n], "med_dg", "med_dv_pub", 3); _, c6, nb6 = rho_block(TB[n], "med_dg", "med_dv_pub", 6)
        rec(f"M2rho|{n}", e, c3, c6, nb3, nb6, None, None)
    AD = pd.read_csv(f"{C03}/resultsD/D_A_tile.csv")
    for n in ("O'E", "N'E"):
        sub = lambda d: AD[(AD.set == n) & (AD["var"] == "pub") & (AD["mode"] == "C1") & (AD["T"] == 5) & (AD.dir == d)]
        for d in ("loss", "gain"):
            e, c3, nb3 = ratio_block(sub(d), "n_lost", "E1", 3); _, c6, nb6 = ratio_block(sub(d), "n_lost", "E1", 6); rec(f"Erho1|{n}|{d}", e, c3, c6, nb3, nb6, None, None)
        s2 = sub("loss"); s2 = s2[s2.n_mov_N2 > 0]; e, c3, nb3 = ratio_block(s2, "n_lost_N2", "E2", 3); _, c6, nb6 = ratio_block(s2, "n_lost_N2", "E2", 6); rec(f"Erho2|{n}", e, c3, c6, nb3, nb6, None, None)
    judged = [k for k, v in r2.items() if isinstance(v, dict) and v.get("verdict_orig") is not None]
    r2["n_judged"] = len(judged); r2["changed"] = [k for k in judged if not r2[k]["unchanged"]]; r2["holds"] = bool(len(r2["changed"]) == 0)
    S["R2"] = r2; pd.DataFrame(rows2).to_csv(f"{OUT}/R2.csv", index=False)
    # ---------------- R3
    cur = []
    for (n, T), g in R3.groupby(["set", "T"]):
        L, G = g[g.dir == "loss"], g[g.dir == "gain"]; L2 = L[L.n_mov_N2 > 0]
        tl = (L.n_lost / L.n_mov.replace(0, np.nan)).dropna()
        cur.append(dict(set=n, T=T, n_mov=int(L.n_mov.sum()), loss_pooled=L.n_lost.sum() / L.n_mov.sum() if L.n_mov.sum() else np.nan,
                        gain_pooled=G.n_lost.sum() / G.n_mov.sum() if G.n_mov.sum() else np.nan, loss_tile_median=float(tl.median()) if len(tl) else np.nan,
                        rho1_loss=L.n_lost.sum() / L.E1.sum() if L.E1.sum() else np.nan, rho1_gain=G.n_lost.sum() / G.E1.sum() if G.E1.sum() else np.nan,
                        rho2=L2.n_lost_N2.sum() / L2.E2.sum() if L2.E2.sum() else np.nan))
    C = pd.DataFrame(cur); C.to_csv(f"{OUT}/R3_curves.csv", index=False)
    chk = {}; Am = {"O": ("O", "C0"), "O'U": ("O'", "C1"), "N'U": ("N'", "C1")}
    for n, (an, md) in Am.items():
        for Tv in (5, 10):
            mine = R3[(R3.set == n) & (R3["T"] == float(Tv))].set_index(["tile", "dir"])[["n_mov", "n_lost", "E1", "E2"]].sort_index()
            ref = A[(A.set == an) & (A["var"] == "pub") & (A["mode"] == md) & (A["T"] == Tv)].set_index(["tile", "dir"])[["n_mov", "n_lost", "E1", "E2"]].sort_index()
            ref = ref.loc[ref.index.intersection(mine.index)]
            chk[f"{n}|T{Tv}"] = bool(len(ref) == len(mine) and len(mine) > 0 and np.allclose(mine.values, ref.loc[mine.index].values, rtol=0, atol=1e-9))
    W = C[C.set.isin(["O", "O'U", "N'U"]) & C["T"].between(3, 12)]
    S["R3"] = dict(R3_0_equal_to_Confirm03=chk, per_set={n: dict(min_rho1_loss=float(g.rho1_loss.min()), min_rho1_gain=float(g.rho1_gain.min()), min_rho2=float(g.rho2.min()),
                                                         all_gt_1=bool((g.rho1_loss > 1).all() and (g.rho1_gain > 1).all() and (g.rho2 > 1).all())) for n, g in W.groupby("set")})
    S["R3"]["holds"] = bool(all(v["all_gt_1"] for v in S["R3"]["per_set"].values()) and len(S["R3"]["per_set"]) == 3)
    # ---------------- R4 (exploratory)
    r4 = {}
    for (n, T), g in R4.groupby(["set", "T"]):
        s = g.sum(numeric_only=True)
        r4[f"{n}|T{int(T)}"] = dict(n_mov=int(s.n_mov), frac_mov_mvs_missing=float(s.n_mov_mvs_missing / s.n_mov), frac_mov_in_buffer=float(s.n_buf / s.n_mov),
                                    frac_losses_in_buffer=float(s.n_lost_buf / s.n_lost) if s.n_lost else np.nan, loss_rate_all=float(s.n_lost / s.n_mov),
                                    loss_rate_outside_buffer=float((s.n_lost - s.n_lost_buf) / (s.n_mov - s.n_buf)) if s.n_mov > s.n_buf else np.nan,
                                    frac_mov_d_lt_1=float(s.n_d1 / s.n_mov), frac_losses_d_lt_1=float(s.n_lost_d1 / s.n_lost) if s.n_lost else np.nan,
                                    loss_rate_outside_d1=float((s.n_lost - s.n_lost_d1) / (s.n_mov - s.n_d1)) if s.n_mov > s.n_d1 else np.nan)
    S["R4"] = r4; pd.DataFrame([dict(key=k, **v) for k, v in r4.items()]).to_csv(f"{OUT}/R4.csv", index=False)
    json.dump(S, open(f"{OUT}/summary.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print("written", OUT, flush=True)
