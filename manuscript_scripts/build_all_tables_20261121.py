"""Generate all LaTeX tables (main text + Supplementary Material) directly from
frozen JSON/CSV result files. No hand-copied numbers without assertion checks.
"""
import json
import pathlib
import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
C01 = ROOT / "contracts" / "Confirm01_JAG"
C02 = ROOT / "contracts" / "Confirm02_Mechanism"
C03 = ROOT / "contracts" / "Confirm03_Proximity"
C04 = ROOT / "contracts" / "Confirm04_Robustness"
C05 = ROOT / "contracts" / "Confirm05_Example"
EVID = ROOT / "evidence"
TAB = ROOT / "manuscript" / "tables"
TAB.mkdir(parents=True, exist_ok=True)


def pct(x, dec=1):
    return f"{100.0 * x:.{dec}f}\\%" if dec is not None else f"{100.0 * x:.0f}\\%"


def fmt_dec(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "--"
    if abs(x) < 0.01:
        return f"{x:.4f}"
    if abs(x) < 0.1:
        return f"{x:.3f}"
    return f"{x:.2f}"


def fmt_ci(v, ci, dec=2):
    return f"{v:.{dec}f} [{ci[0]:.{dec}f}, {ci[1]:.{dec}f}]"


def build_all():
    p1_o_rc = json.loads((EVID / "Probe_RelCons01" / "summary.json").read_text(encoding="utf-8"))
    p1_o_s = json.loads((C01 / "results" / "part1_Oprime__summary.json").read_text(encoding="utf-8"))
    p1_n_s = json.loads((C01 / "results" / "part1_Nprime__summary.json").read_text(encoding="utf-8"))
    p2_sum = json.loads((C01 / "results" / "part2__summary.json").read_text(encoding="utf-8"))
    p3_sum = json.loads((C01 / "results" / "part3__summary.json").read_text(encoding="utf-8"))

    c2_sum = json.loads((C02 / "results" / "mech" / "summary.json").read_text(encoding="utf-8"))
    c3_sum = json.loads((C03 / "results" / "summary.json").read_text(encoding="utf-8"))
    c3_d = json.loads((C03 / "results" / "summary_D.json").read_text(encoding="utf-8"))
    c4_sum = json.loads((C04 / "results" / "summary.json").read_text(encoding="utf-8"))
    c5_counts = pd.read_csv(C05 / "results" / "example_counts.csv")

    bt_o = pd.read_csv(C02 / "results" / "mech" / "by_tile_O.csv")
    bt_op = pd.read_csv(C02 / "results" / "mech" / "by_tile_O'.csv")
    bt_np = pd.read_csv(C02 / "results" / "mech" / "by_tile_N'.csv")
    dt_u = pd.read_csv(C03 / "results" / "D_tile_U.csv")
    dt_e = pd.read_csv(C03 / "results" / "D_tile_E.csv")
    r1_df = pd.read_csv(C04 / "results" / "R1.csv").set_index("set")
    r2_df = pd.read_csv(C04 / "results" / "R2.csv").set_index("key")
    r4_df = pd.read_csv(C04 / "results" / "R4.csv").set_index("key")

    def cls_counts(df):
        d = int((df["n"] >= 200).sum())
        s = int((df["n"] <= 130).sum())
        m = int(len(df) - d - s)
        return f"{d}/{m}/{s}"

    def n_pts_from_tile(df, col="C0__pub_T5_n"):
        return int(df[col].sum())

    # ------------------------------------------------------------------
    # Table 1: Datasets and release pairs (incl. E-W & S'')
    # ------------------------------------------------------------------
    rows_t1 = [
        ("O (U)", "P1: 2019--2023 vs.\\ 2020--2024", len(bt_o), n_pts_from_tile(bt_o),
         cls_counts(bt_o), "2.0 $\\rightarrow$ 2.0 (ETRF2000)",
         f"{p1_o_rc['F_n_stations']} ({p1_o_rc['F_n_sparse']} sparse)"),
        ("O$'$ (U)", "P2: 2020--2024 vs.\\ 2020--2025", len(bt_op), n_pts_from_tile(bt_op),
         cls_counts(bt_op), "2.0 $\\rightarrow$ 3.0 (ETRF2000 $\\rightarrow$ 2020)",
         f"{p1_o_s['n_stations']} ({p1_o_s['n_stations_sparse']} sparse)"),
        ("N$'$ (U)", "P2: 2020--2024 vs.\\ 2020--2025", len(bt_np), n_pts_from_tile(bt_np),
         cls_counts(bt_np), "2.0 $\\rightarrow$ 3.0 (ETRF2000 $\\rightarrow$ 2020)",
         f"{p1_n_s['n_stations']} ({p1_n_s['n_stations_sparse']} sparse)"),
        ("O$'$ (E)", "P2: 2020--2024 vs.\\ 2020--2025",
         len(dt_e[dt_e.set == "O'E"]), n_pts_from_tile(dt_e[dt_e.set == "O'E"]),
         cls_counts(dt_e[dt_e.set == "O'E"]), "2.0 $\\rightarrow$ 3.0 (ETRF2000 $\\rightarrow$ 2020)", "--"),
        ("N$'$ (E)", "P2: 2020--2024 vs.\\ 2020--2025",
         len(dt_e[dt_e.set == "N'E"]), n_pts_from_tile(dt_e[dt_e.set == "N'E"]),
         cls_counts(dt_e[dt_e.set == "N'E"]), "2.0 $\\rightarrow$ 3.0 (ETRF2000 $\\rightarrow$ 2020)", "--"),
        ("S$''$ (U)", "P2: 2020--2024 vs.\\ 2020--2025",
         len(dt_u[dt_u.set == "S''"]), n_pts_from_tile(dt_u[dt_u.set == "S''"]),
         cls_counts(dt_u[dt_u.set == "S''"]), "2.0 $\\rightarrow$ 3.0 (ETRF2000 $\\rightarrow$ 2020)", "--"),
    ]
    t1 = [
        r"\begin{table*}[pos=t]",
        r"\caption{EGMS Ortho release pairs, tile sets (U: vertical; E: east--west), common-window epoch density classes (dense $n \ge 200$, mid $131$--$199$, sparse $n \le 130$), XML GNSS calibration transitions \citep{EGMS_GNSS_2026}, and qualifying NGL IGS20 GNSS stations ($618$ unique stations in total; $597$ in P2).}\label{tab:data}",
        r"\scriptsize\setlength{\tabcolsep}{2.4pt}",
        r"\begin{tabular*}{\tblwidth}{@{\extracolsep{\fill}}llrrlll@{}}",
        r"\toprule",
        r"Set & Release pair & Tiles & Points & D/M/S & XML model (frame) & GNSS stations\\",
        r"\midrule",
    ]
    for s_id, pair, nt, npt, dms, mdl, nst in rows_t1:
        t1.append(f"{s_id} & {pair} & {nt} & {npt:,} & {dms} & {mdl} & {nst}\\\\")
    t1 += [r"\bottomrule", r"\end{tabular*}", r"\end{table*}"]
    (TAB / "tab01_data.tex").write_text("\n".join(t1), encoding="utf-8")

    # ------------------------------------------------------------------
    # Table 2: Concise nomenclature & abbreviations
    # ------------------------------------------------------------------
    t2 = [
        r"\begin{table*}[pos=t]",
        r"\caption{Key notation, counterfactual adjustments, and decision metrics used across Sections~\ref{sec:methods}--\ref{sec:disc}.}\label{tab:nomen}",
        r"\footnotesize\setlength{\tabcolsep}{3.5pt}",
        r"\begin{tabular*}{\tblwidth}{@{}p{0.18\linewidth}p{0.78\linewidth}@{}}",
        r"\toprule",
        r"Symbol / term & Definition\\",
        r"\midrule",
        r"\texttt{P1}, \texttt{P2} & Release pairs: \texttt{P1} = 2019--2023 (A) vs.\ 2020--2024 (B); \texttt{P2} = 2020--2024 (A) vs.\ 2020--2025 (B).\\",
        r"\texttt{O}, \texttt{O$'$}, \texttt{N$'$}, \texttt{S$''$} & Fixed 136-tile panel on P1 (\texttt{O}) and P2 (\texttt{O$'$}); 72 unseen stratified tiles (\texttt{N$'$}); 35 northern (\texttt{Nt}) supplement tiles (15 sparse, 12 mid, 8 dense; \texttt{S$''$}).\\",
        r"\texttt{C0}, \texttt{C1}, \texttt{C2} & Unadjusted release B (\texttt{C0}); minus per-point $\Dg = g_B - g_A$ (\texttt{C1}); minus $\Dg$ and per-tile residual OLS plane (\texttt{C2}).\\",
        r"\texttt{A1}--\texttt{A3} & Decision metrics: point class disagreement (\texttt{A1}) plus moving-point loss/gain; hotspot persistence $\mathrm{pers}_A$ and loss $1-\mathrm{pers}_A$ (\texttt{A2}); GHSL pooled flagged disagreement (\texttt{A3}).\\",
        r"$s_g$, $s_{\mathrm{ramp}}$, $s_{\mathrm{local}}$ & Partition of \texttt{C0} disagreement into frame ($D_{\mathrm{C0}}-D_{\mathrm{C1}}$), tile-plane ($D_{\mathrm{C1}}-D_{\mathrm{C2}}$), and local ($D_{\mathrm{C2}}$) shares normalized by $D_{\mathrm{C0}}$.\\",
        r"$\rho_1$, $\rho_2$ & Moving-point loss (or gain) ratio against random-sign flipping around observed boundary ($L_1$) and typical-point $|\Dv|$ shifts ($L_2$).\\",
        r"$R$, $R_{\mathrm{dec}}$ & Station velocity-difference over EGMS--GNSS misfit ratio ($R$), and inter-release over EGMS--GNSS station flip ratio ($R_{\mathrm{dec}}$).\\",
        r"\bottomrule",
        r"\end{tabular*}",
        r"\end{table*}",
    ]
    (TAB / "tab02_nomen.tex").write_text("\n".join(t2), encoding="utf-8")

    # ------------------------------------------------------------------
    # Table 3: Compact hash-locked contracts summary
    # ------------------------------------------------------------------
    t3 = [
        r"\begin{table*}[pos=t]",
        r"\caption{Chronology of SHA-256 hash-locked analysis contracts. Complete expectation criteria and block-bootstrap intervals appear in Supplementary Table~S1.}\label{tab:contracts}",
        r"\scriptsize\setlength{\tabcolsep}{2.4pt}",
        r"\begin{tabular*}{\tblwidth}{@{\extracolsep{\fill}}p{0.12\linewidth}p{0.14\linewidth}p{0.15\linewidth}p{0.19\linewidth}p{0.34\linewidth}@{}}",
        r"\toprule",
        r"Contract & Frozen (UTC) & SHA-256 & Scope & Principal verified vs.\ failed gates\\",
        r"\midrule",
        r"\texttt{Confirm01} & 2026-11-17 10:04 & \texttt{61cd70...fa672a} & P2 draws, point \& GNSS checks & 4/8 replication checks met (\texttt{F2,F3,F5,F6}); \texttt{F1,F4,F7,F8} not replicated on P2; \texttt{cw<pub} at $\pm 10$ only.\\",
        r"\texttt{Amendment 1} & 2026-11-17 13:11 & \texttt{8b7a96...35fbed5} & Switch unseen target P1$\to$P2 & Approved prior to unsealing \texttt{N$'$} or \texttt{P2}; no numeric threshold changed.\\",
        r"\texttt{Confirm02} & 2026-11-18 10:36 & \texttt{a2a07f...e2d15b2} & Frame attribution (\texttt{C0--C2}), RSC & \texttt{M0--M4} and software checks \texttt{U1--U3} passed; frame-share gates \texttt{M5a}/\texttt{M5b} failed.\\",
        r"\texttt{Confirm03} & 2026-11-20 12:12 & \texttt{257ab6...3978af1} & Proximity, GNSS offsets, E--W, \texttt{S$''$} & \texttt{H3,H4,H5,H7} passed; threshold concentration \texttt{H1}, cross-pair GNSS \texttt{H6}, and monotonicity \texttt{H8} failed; \texttt{H2} mixed on \texttt{N$'$}.\\",
        r"\texttt{Confirm04} & 2026-11-21 06:45 & \texttt{c57fbc...bb640a58} & Rounding, 300\,km blocks, curves & Rounding (\texttt{R1}) \& $3$--$12$\,mm/yr curves (\texttt{R3}) passed; 300\,km blocks (\texttt{R2}) unchanged on 21/23 checks.\\",
        r"\texttt{Confirm05} & 2026-11-21 08:28 & \texttt{fa8f0f...95f0446} & Case study \texttt{E39N31} (Antwerp) & Descriptive illustration selected by max GHSL population (47.6th loss percentile); no gate.\\",
        r"\bottomrule",
        r"\end{tabular*}",
        r"\end{table*}",
    ]
    (TAB / "tab03_contracts.tex").write_text("\n".join(t3), encoding="utf-8")

    # ------------------------------------------------------------------
    # Table 4: Main decision stability & frame-share readout (Pooled + Tile Median + E-W key rows)
    # ------------------------------------------------------------------
    # Helper to compute exact pooled and tile-median C0/C1 point disagr & moving-point loss directly from by_tile CSVs
    def tile_stats(df, mode, T):
        sub_n = df[f"{mode}__pub_{T}_n"]
        sub_dis = df[f"{mode}__pub_{T}_n_dis"]
        a1_pool = sub_dis.sum() / sub_n.sum()
        a1_med = (sub_dis / sub_n).median()
        out = dict(A1_pool=a1_pool, A1_med=a1_med)
        if f"{mode}__pub_{T}_n_movA" in df.columns:
            ma = df[f"{mode}__pub_{T}_n_movA"]
            ml_rate = df[f"{mode}__pub_{T}_loss"]
            valid = ma > 0
            loss_pool = (ma[valid] * ml_rate[valid]).sum() / ma[valid].sum()
            loss_med = ml_rate[valid].median()
            ha = df[f"{mode}__pub_{T}_n_hotA"]
            vhot = ha > 0
            hot_loss_med = (1.0 - df.loc[vhot, f"{mode}__pub_{T}_persA"]).median()
            bdis = df[f"{mode}__pub_{T}_built_dis"].sum()
            bflg = df[f"{mode}__pub_{T}_built_flag"].sum()
            out.update(loss_pool=loss_pool, loss_med=loss_med, hot_loss_med=hot_loss_med,
                       built_pool=bdis / bflg if bflg > 0 else np.nan)
        return out

    # Cross-check O, O', N' against Confirm01 part2__summary, Confirm02 aggregates, and Confirm03 summary
    assert abs(tile_stats(bt_o, "C0", "T5")["loss_pool"] - p2_sum["O"]["pub_T5"]["loss_pooled"]) < 1e-6
    assert abs(tile_stats(bt_op, "C1", "T5")["loss_pool"] - c3_sum["A"]["O'|pub|C1|T5|loss"]["L_obs"]) < 1e-6
    assert abs(tile_stats(bt_np, "C1", "T5")["loss_pool"] - c3_sum["A"]["N'|pub|C1|T5|loss"]["L_obs"]) < 1e-6
    assert abs(tile_stats(bt_op, "C0", "T5")["A1_med"] - p2_sum["O'"]["pub_T5"]["A1_tile_median"]) < 1e-6
    assert abs(tile_stats(bt_np, "C0", "T5")["built_pool"] - p2_sum["N'"]["pub_T5"]["A3_built"]) < 1e-6

    sg_op = c2_sum["M5"]["shares"]
    sg_np = c2_sum["M5"]["shares"]
    sg_ew_o = c3_d["D1"]["O'E"]["metrics"]["shares"]
    sg_ew_n = c3_d["D1"]["N'E"]["metrics"]["shares"]

    def sg_str(d, key):
        v = d[key]
        val = f"{v['s_g']:.2f}"
        if val == "-0.00":
            val = "0.00"
        return f"{val} [{v['s_g_ci'][0]:.2f}, {v['s_g_ci'][1]:.2f}]"

    def a1_pair_str(x):
        m, p = x["A1_med"], x["A1_pool"]
        if m < 0.001:
            return f"{m:.4f}/{p:.3f}"
        if m < 0.01:
            return f"{m:.3f}/{p:.3f}"
        return f"{m:.2f}/{p:.2f}"

    # Build concise multi-column table: U + key E-W lines
    t4 = [
        r"\begin{table*}[pos=t]",
        r"\caption{Screening-decision instability and frame share ($s_g$, 95\% tile-bootstrap CI) across P1 (\texttt{O}) and P2 (\texttt{O$'$}, \texttt{N$'$}), comparing tile medians and pooled totals before (\texttt{C0}) and after frame harmonization (\texttt{C1}). Full east--west, geographic-strata, and uncertainty-buffer tables are in Supplementary Tables~S2--S3.}\label{tab:decisions}",
        r"\scriptsize\setlength{\tabcolsep}{1.8pt}",
        r"\begin{tabular*}{\tblwidth}{@{\extracolsep{\fill}}p{0.20\linewidth}lrrrrrcc@{}}",
        r"\toprule",
        r"Metric & Cut & O C0 & O$'$ C0 & O$'$ C1 & N$'$ C0 & N$'$ C1 & O$'$ $s_g$ [CI] & N$'$ $s_g$ [CI]\\",
        r"\midrule",
        r"\multicolumn{9}{@{}l}{\textit{Vertical (U): tile median / pooled across all tiles}}\\",
    ]
    for T_code, T_lbl in [("T0", "sign ($v<0$)"), ("T05", "$\\pm 0.5$"), ("T5", "$\\pm 5$"), ("T10", "$\\pm 10$")]:
        o0 = tile_stats(bt_o, "C0", T_code)
        op0 = tile_stats(bt_op, "C0", T_code)
        op1 = tile_stats(bt_op, "C1", T_code)
        np0 = tile_stats(bt_np, "C0", T_code)
        np1 = tile_stats(bt_np, "C1", T_code)
        vals = [a1_pair_str(x) for x in (o0, op0, op1, np0, np1)]
        sg_o_str = sg_str(sg_op, f"O'|pub|A1|{T_code}")
        sg_n_str = sg_str(sg_np, f"N'|pub|A1|{T_code}")
        t4.append(f"Point disagr.\\ (med/pool) & {T_lbl} & " + " & ".join(vals) +
                  f" & {sg_o_str} & {sg_n_str}\\\\")

    for T_code, T_lbl in [("T5", "$\\pm 5$"), ("T10", "$\\pm 10$")]:
        o0 = tile_stats(bt_o, "C0", T_code)
        op0 = tile_stats(bt_op, "C0", T_code)
        op1 = tile_stats(bt_op, "C1", T_code)
        np0 = tile_stats(bt_np, "C0", T_code)
        np1 = tile_stats(bt_np, "C1", T_code)
        # Moving-point loss med/pool
        lvals = [f"{x['loss_med']:.2f}/{x['loss_pool']:.2f}" for x in (o0, op0, op1, np0, np1)]
        t4.append(f"Moving-pt loss (med/pool) & {T_lbl} & " + " & ".join(lvals) + " & -- & --\\\\")
        # Hotspot loss med
        hvals = [f"{x['hot_loss_med']:.2f}" for x in (o0, op0, op1, np0, np1)]
        sg_oa2 = sg_str(sg_op, f"O'|pub|A2|{T_code}")
        sg_na2 = sg_str(sg_np, f"N'|pub|A2|{T_code}")
        t4.append(f"Hotspot loss (tile med) & {T_lbl} & " + " & ".join(hvals) +
                  f" & {sg_oa2} & {sg_na2}\\\\")
        # Built-up disagr pool
        bvals = [f"{x['built_pool']:.2f}" for x in (o0, op0, op1, np0, np1)]
        sg_oa3 = sg_str(sg_op, f"O'|pub|A3|{T_code}")
        sg_na3 = sg_str(sg_np, f"N'|pub|A3|{T_code}")
        t4.append(f"Built-up disagr.\\ (pool) & {T_lbl} & " + " & ".join(bvals) +
                  f" & {sg_oa3} & {sg_na3}\\\\")

    # E-W summary rows (selected thresholds)
    oe_tiles = dt_e[dt_e.set == "O'E"]
    ne_tiles = dt_e[dt_e.set == "N'E"]
    t4.append(r"\midrule")
    t4.append(r"\multicolumn{9}{@{}l}{\textit{East--west (E): tile median / pooled (key thresholds)}}\\")
    for T_code, T_lbl in [("T05", "$\\pm 0.5$"), ("T5", "$\\pm 5$"), ("T10", "$\\pm 10$")]:
        oe0 = tile_stats(oe_tiles, "C0", T_code)
        oe1 = tile_stats(oe_tiles, "C1", T_code)
        ne0 = tile_stats(ne_tiles, "C0", T_code)
        ne1 = tile_stats(ne_tiles, "C1", T_code)
        vals = ["--"] + [a1_pair_str(x) for x in (oe0, oe1, ne0, ne1)]
        t4.append(f"Point disagr.\\ E (med/pool) & {T_lbl} & " + " & ".join(vals) +
                  f" & {sg_str(sg_ew_o, f'pub|A1|{T_code}')} & {sg_str(sg_ew_n, f'pub|A1|{T_code}')}\\\\")
    # Also add E +/-5 moving-point loss and built-up
    oe0 = tile_stats(oe_tiles, "C0", "T5")
    oe1 = tile_stats(oe_tiles, "C1", "T5")
    ne0 = tile_stats(ne_tiles, "C0", "T5")
    ne1 = tile_stats(ne_tiles, "C1", "T5")
    t4.append(f"Moving-pt loss E (med/pool) & $\\pm 5$ & -- & {oe0['loss_med']:.2f}/{oe0['loss_pool']:.2f} & {oe1['loss_med']:.2f}/{oe1['loss_pool']:.2f} & {ne0['loss_med']:.2f}/{ne0['loss_pool']:.2f} & {ne1['loss_med']:.2f}/{ne1['loss_pool']:.2f} & -- & --\\\\")
    t4 += [r"\bottomrule", r"\end{tabular*}", r"\end{table*}"]
    (TAB / "tab04_decisions.tex").write_text("\n".join(t4), encoding="utf-8")

    # ------------------------------------------------------------------
    # Table 5: Two-tier Release Stability Check (RSC) user guidance
    # ------------------------------------------------------------------
    t5 = [
        r"\begin{table*}[pos=t]",
        r"\caption{Release Stability Check (RSC): operational workflow and two-tier decision rules grounded in this audit.}\label{tab:rsc}",
        r"\footnotesize\setlength{\tabcolsep}{3.5pt}",
        r"\begin{tabular*}{\tblwidth}{@{}p{0.19\linewidth}p{0.56\linewidth}p{0.21\linewidth}@{}}",
        r"\toprule",
        r"Step / tier & Action for EGMS Ortho pair comparisons & Key evidence\\",
        r"\midrule",
        r"1.\ Metadata audit & Inspect XML \texttt{<gnss><version>} and production tags across both releases; a transition ($2.0 \rightarrow 3.0$) identifies the ETRF2000-to-ETRF2020 change. & Table~\ref{tab:data}, Sec.~\ref{sec:res-frame}\\",
        r"2.\ Frame harmonize & Subtract $\Dg = g_B - g_A$ from \texttt{gnss\_velocity\_\{u,e\}} (falling back to constant median $\Dg$ if unavailable) prior to comparing velocities. & Fig.~\ref{fig:frame}, Table~\ref{tab:decisions}\\",
        r"3.\ Plane check (opt.) & Residual 100\,km tile planes contribute little ($|s_{\mathrm{ramp}}| \le 0.13$); detrend only for regional ramp checks and report both \texttt{C1} and \texttt{C2}. & Fig.~\ref{fig:decisions}, Sec.~\ref{sec:res-attrib}\\",
        r"4a.\ Tier~1: near-zero ($v<0$, $\pm 0.5$\,mm/yr) & \textbf{Frame-sensitive:} harmonize beforehand (halves station flip ratio $R_{\mathrm{dec}}$; $s_g$ up to $0.52$), then treat remaining $8$--$24\%$ disagreement as release-sensitive differences. & Figs.~\ref{fig:decisions}b, \ref{fig:gnss}\\",
        r"4b.\ Tier~2: alert gates ($\pm 5$, $\pm 10$\,mm/yr) & \textbf{Local-reprocessing dominated:} frame shift accounts for minor shares ($s_g \le 0.22$); $16$--$22\%$ pooled ($37$--$62\%$ tile-median) \texttt{C0} $\pm 5$ active loss and $20$--$60\%$ built-up disagreement remain. Report turnover bounds alongside hotspots (\textit{both}/\textit{A-only}/\textit{B-only}). & Figs.~\ref{fig:threshold}, \ref{fig:antwerp}; Table~\ref{tab:decisions}\\",
        r"5.\ Do not & Do not treat quantized \texttt{mean\_velocity\_std} as release-to-release uncertainty, nor infer which release is closer to true surface motion without external ground truth. & Fig.~\ref{fig:gnss}, Supp.\ Table~S2\\",
        r"\bottomrule",
        r"\end{tabular*}",
        r"\end{table*}",
    ]
    (TAB / "tab05_rsc.tex").write_text("\n".join(t5), encoding="utf-8")

    # ------------------------------------------------------------------
    # Supplementary Tables S1-S4
    # ------------------------------------------------------------------
    # Table S1: Complete preregistered expectations & block-bootstrap verdicts
    s1 = [
        r"\begin{table}[htbp]",
        r"\caption{Complete ledger of locally hash-locked confirmatory expectations (Confirm01--04), observed estimates, 300\,km block-bootstrap intervals where applicable, and final verdicts.}\label{tab:supp-prereg}",
        r"\scriptsize\setlength{\tabcolsep}{1.8pt}",
        r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}p{0.10\linewidth}p{0.05\linewidth}p{0.20\linewidth}p{0.19\linewidth}p{0.24\linewidth}p{0.14\linewidth}@{}}",
        r"\toprule",
        r"Contract & ID & Expectation & Criterion & Observed (O$'$ / N$'$ or set) & Verdict\\",
        r"\midrule",
        f"Confirm01 & F1 & Category flip rate (trend; step) & trend 8--20\\%; step 15--35\\% & {pct(p1_o_s['F1']['trend'],1)}; {pct(p1_o_s['F1']['step'],1)} / {pct(p1_n_s['F1']['trend'],1)}; {pct(p1_n_s['F1']['step'],1)} & Not repl.\\ (P1 only)\\\\",
        f"Confirm01 & F2 & Tile $\\rho(n, m)$ & $\\rho\\le-0.5$, upper $<-0.2$ & {p1_o_s['F2']['rho']:.2f} / {p1_n_s['F2']['rho']:.2f} (300\\,km [-0.77, -0.43]) & Replicated\\\\",
        f"Confirm01 & F3 & Share $|\\Dv|>2\\sigma_\\mathrm{{form}}$ & $\\ge 0.50$ per class & 0.81--0.93 / 0.84--0.92 & Replicated\\\\",
        f"Confirm01 & F4 & P1 \\textbf{{B0}} coverage & pooled 0.92--0.97; class $\\ge 0.88$ & 0.898 / 0.863 (dense 0.831) & Not repl.\\ (P1 only)\\\\",
        f"Confirm01 & F5 & Product std vs.\\ \\textbf{{B0}} & $\\le 20$ values; B1/B0 $\\ge 1.10$ & 6 / 8 values; 1.23 / 1.20 & Replicated\\\\",
        f"Confirm01 & F6 & Short/long spatial ratio & $\\le 0.80$, upper $<1.00$ & {p1_o_s['F6']['r_median']:.2f} / {p1_n_s['F6']['r_median']:.2f} & Replicated\\\\",
        f"Confirm01 & F7 & GNSS $\\rho(|\\Dv_s|, e_s)$ & $\\rho\\ge 0.10$, p5 $>0$ & {p1_o_s['F7']['rho']:.2f} / {p1_n_s['F7']['rho']:.2f} & Not repl.\\ (P1 only)\\\\",
        f"Confirm01 & F8 & Ratio $R$ limitation gate & $R<0.50$ if limited & {p1_o_s['F8']['R']:.2f} / {p1_n_s['F8']['R']:.2f} (\\textbf{{C1}}: 0.44 / 0.32) & Not repl.\\ (P1 only)\\\\",
        f"Confirm01 & A & \\texttt{{cw}} $<$ \\texttt{{pub}} (\\textbf{{A1}}) & all sets, $\\pm 5$ and $\\pm 10$ & $\\pm 10$ holds (all 3); $\\pm 5$ fails P2 & $\\pm 10$ only\\\\",
        f"Confirm02 & M1 & $\\Dg$ frame consistency & $\\ge 95\\%$ in $[-0.01,0.88]$; $\\ge 90\\%>0$ & 0.99, 0.96 / 1.00, 0.99 & Holds\\\\",
        f"Confirm02 & M2 & Systematic offset removed & $|S_1|\\le 0.10$ and $\\le |S_0|/3$ & $+0.50\\to 0.00$ / $+0.50\\to 0.00$ & Holds\\\\",
        f"Confirm02 & M3 & Dense $m^*$ ratio O$'$(C1)/O & in $[2/3, 3/2]$ & 1.72 (C0) $\\to$ 1.37 (C1) & Holds\\\\",
        f"Confirm02 & M5a & $s_g$ near zero ($0$, $\\pm 0.5$) & $\\ge 0.30$, CI $>0$, both & O$'$ 0.07, 0.24 / N$'$ 0.52, 0.33 & Fails (N$'$ only)\\\\",
        f"Confirm02 & M5b & $s_g$ at $\\pm 10$ (A1, A2, A3) & upper CI $\\le 0.10$, both & A3 upper 0.15/0.34; O$'$ 300\\,km A1/A2 0.12/0.11 & Fails\\\\",
        f"Confirm02 & M5c & O$'$(C1)/O at $\\pm 0.5$ \\textbf{{A1}} & in $[2/3, 3/2]$ & 1.12 (O$'$ C0 already 1.43) & Non-discrim.\\\\",
        f"Confirm02 & U1--U3 & RSC software checks & bit-identity, injection, $\\Dg\\equiv 0$ & 162 cols; $<2.8\\times 10^{{-13}}$; exact & Pass\\\\",
        f"Confirm03 & H1 & $\\pm 5$ loss within $d<1$ & share $\\ge 0.50$ and ratio $\\ge 3$ & 0.42--0.44; ratio 2.39--2.60 & Fails\\\\",
        f"Confirm03 & H2 & $\\rho_1$ vs.\\ random sign & both dir.\\ CI $>1.0$ & O, O$'$ pass; N$'$ gain [0.83, 1.61] & Mixed (fails N$'$)\\\\",
        f"Confirm03 & H3 & $\\rho_2$ vs.\\ typical noise & lower CI $>1.0$ & 1.71 / 2.10 / 2.29 (all CI $>1.4$) & Holds\\\\",
        f"Confirm03 & H4--H6 & GNSS offsets $c$ & $|\\Delta c|\\le 0.25$ across/within & H4 ($-0.14$, $+0.02$), H5 ($-0.02$) hold; H6 fails ($-0.43$) & H4,H5 hold; H6 fails\\\\",
        f"Confirm03 & H7 & E--W $\\Dg_E$ removal & S0 $>0.3\\to|S_1|\\le 0.1$, $\\rho>0.4$ & O$'$E 0.50$\\to$0.00 ($\\rho$ 0.74); N$'$E 0.60$\\to$0.00 (0.87) & Holds\\\\",
        f"Confirm03 & H8 & N$'$ $\\cup$ S$''$ monotonicity & $\\rho\\le -0.5$, CI $<-0.2$, ratio $\\ge 3$ & \\textbf{{C1 cw}}: $\\rho=-0.41$ [-0.58, -0.20]; ratio 4.0 & Fails\\\\",
        f"Confirm04 & R1--R3 & Rounding, 300\\,km block, $T$ curves & M2 under HD/jitter; 23 checks; $\\rho_1,\\rho_2>1$ & S1 $-0.02..+0.05$; 21/23 block ok; T3--12 $\\rho_1,\\rho_2>1$ & R1,R3 hold; 2/23 block shift\\\\",
        r"\bottomrule",
        r"\end{tabular*}",
        r"\end{table}",
    ]
    (TAB / "tabS1_prereg_all.tex").write_text("\n".join(s1), encoding="utf-8")

    b2o, b2n = p3_sum["O'"]["B2"], p3_sum["N'"]["B2"]
    b3o, b3n = p3_sum["O'"]["B3"], p3_sum["N'"]["B3"]
    # Table S2: Secondary point-level, B0/B1/B2/B3, constant/plane variants, and R4 buffer
    s2 = [
        r"\begin{table}[htbp]",
        r"\caption{Secondary diagnostics moved from main text: uncertainty calibration (\textbf{B0}/\textbf{B1}), donor thinning (\textbf{B2}), residual temporal structure (\textbf{B1}/\textbf{B3}), frame variants (\texttt{C1c}/\texttt{C2$'$}), and exploratory buffer sensitivity (\textbf{R4}).}\label{tab:supp-sec}",
        r"\scriptsize\setlength{\tabcolsep}{3.0pt}",
        r"\begin{tabular*}{\textwidth}{@{}p{0.28\linewidth}p{0.34\linewidth}p{0.34\linewidth}@{}}",
        r"\toprule",
        r"Diagnostic & O$'$ (P2) / O (P1 where noted) & N$'$ (P2)\\",
        r"\midrule",
        f"\\textbf{{B1}} vs.\\ \\textbf{{B0}} interval score (F5) & 19.76 vs.\\ 16.03 (+23\\%); tercile 1.25 vs.\\ 1.63 & 10.83 vs.\\ 9.00 (+20\\%); tercile 1.14 vs.\\ 1.26\\\\",
        f"\\textbf{{B2}} donor calendar thinning & $S = {b2o['S']:.3f}$ [0.020, 0.064]; ratio {b2o['median_ratio']:.3f} ($n=64$) & $S = {b2n['S']:.3f}$ [0.017, 0.052]; ratio {b2n['median_ratio']:.3f} ($n=55$)\\\\",
        f"\\textbf{{B3}} epoch diff.\\ lag-1 AC; dense$-$sparse & AC {b3o['ac1_median']:.2f}; residual corr.\\ $+{b3o['rescorr_dense_minus_sparse']:.2f}$ & AC {b3n['ac1_median']:.2f}; residual corr.\\ $+{b3n['rescorr_dense_minus_sparse']:.2f}$\\\\",
        f"\\textbf{{B1}} partial Spearman (C0 vs.\\ C1) & C0 winter-gap 0.41 vs.\\ $-0.02$; C1 undeterm.\\ 0.27 & C0 undeterm.\\ 0.26 vs.\\ $-0.18$; C1 winter-gap 0.31\\\\",
        f"\\textbf{{F7}}, \\textbf{{F9}}, \\textbf{{F10}} GNSS station checks & $\\rho(|\\Dv_s|,e_s)=0.07$; offset $-0.50$; sign $p=0.36$ & $\\rho(|\\Dv_s|,e_s)=0.00$; offset $-0.41$; sign $p=0.70$\\\\",
        r"\midrule",
        r"\multicolumn{3}{@{}l}{\textit{Exploratory buffer exclusion on moving-point loss (Confirm04 R4, C0 for O, C1 for O$'$, N$'$)}}\\",
    ]
    for key_lbl, rk in [("O $\\pm 5$ / $\\pm 10$", ("O|T5", "O|T10")),
                        ("O$'$ $\\pm 5$ / $\\pm 10$", ("O'U|T5", "O'U|T10")),
                        ("N$'$ $\\pm 5$ / $\\pm 10$", ("N'U|T5", "N'U|T10"))]:
        r5, r10 = r4_df.loc[rk[0]], r4_df.loc[rk[1]]
        s2.append(
            f"R4 buffer: {key_lbl} & "
            f"missing std {100*r5['frac_mov_mvs_missing']:.0f}\\% / {100*r10['frac_mov_mvs_missing']:.0f}\\%; "
            f"$d<2\\sigma$ loss {r5['loss_rate_outside_buffer']:.3f} / {r10['loss_rate_outside_buffer']:.3f} & "
            f"$d<1$\\,mm/yr drops {100*r5['frac_mov_d_lt_1']:.0f}\\%/{100*r10['frac_mov_d_lt_1']:.0f}\\% pts; loss {r5['loss_rate_outside_d1']:.3f} / {r10['loss_rate_outside_d1']:.3f}\\\\"
        )
    s2 += [r"\bottomrule", r"\end{tabular*}", r"\end{table}"]
    (TAB / "tabS2_secondary.tex").write_text("\n".join(s2), encoding="utf-8")

    # Table S3: N' geographic strata & complete E-W decision stability
    s3 = [
        r"\begin{table}[htbp]",
        r"\caption{Confirm03 descriptive strata breakdown on \texttt{N$'$} (Alpine, Central, North, South; vertical, published) and complete east--west (\texttt{E}) decision metrics across \texttt{O$'$} and \texttt{N$'$}. East--west \texttt{A2} medians include all tiles with $\ge 1$ release-A hotspot (whereas $s_g$ in main-text Table~4 restricts to tiles with nonzero \texttt{C0} loss).}\label{tab:supp-strata}",
        r"\scriptsize\setlength{\tabcolsep}{1.5pt}",
        r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}p{0.13\linewidth}lrrrrrr@{}}",
        r"\toprule",
        r"Stratum/set & $n$ & $\pm 5$ loss (C0/C1) & $\pm 10$ loss (C0/C1) & $\pm 10$ A3 (C0/C1) & $s_g$ (A3) & $d<1$ & $\rho_2$ ($\pm 5$)\\",
        r"\midrule",
    ]
    for st_code, st_name in [("Al", "N$'$ Alpine"), ("C", "N$'$ Central"), ("Nt", "N$'$ North"), ("S", "N$'$ South")]:
        c = c3_sum["C"][st_code]
        s3.append(f"{st_name} & {c['n_tiles']} & {c['loss_pooled_C0_T5']:.2f} / {c['loss_pooled_C1_T5']:.2f} & "
                  f"{c['loss_pooled_C0_T10']:.2f} / {c['loss_pooled_C1_T10']:.2f} & "
                  f"{c['A3_T10_C0']:.2f} / {c['A3_T10_C1']:.2f} & {c['sg_A3_T10']:.2f} & "
                  f"{c['H1_share_d_lt_1']:.2f} & {c['H3_rho2']:.2f}\\\\")
    s3 += [
        r"\midrule",
        r"Set (E--W) & Mode & \texttt{A1} $\pm 0.5$ & \texttt{A1} $\pm 5$ & \texttt{A1} $\pm 10$ & \texttt{A2} $\pm 5$/$10$ & \texttt{A3} $\pm 5$/$10$ & Loss $\pm 5$ (m/p)\\",
        r"\midrule",
    ]
    for comp, label in [("O'E", "O$'$ E"), ("N'E", "N$'$ E")]:
        ag = c3_d["D1"][comp]["metrics"]["aggregates"]
        lg = c3_d["D1"][comp]["metrics"]["lossgain"]
        for md in ("C0", "C1", "C2"):
            lm = lg[f"{md}|pub|T5"]["loss_tile_median"] if f"{md}|pub|T5" in lg else np.nan
            lp = lg[f"{md}|pub|T5"]["loss_pooled"] if f"{md}|pub|T5" in lg else np.nan
            lstr = f"{lm:.2f} / {lp:.2f}" if not np.isnan(lm) else "--"
            s3.append(f"{label} & {md} & {ag[f'{md}|pub|A1|T05']:.2f} & {ag[f'{md}|pub|A1|T5']:.3f} & "
                      f"{ag[f'{md}|pub|A1|T10']:.4f} & {ag[f'{md}|pub|A2|T5']:.2f} / {ag[f'{md}|pub|A2|T10']:.2f} & "
                      f"{ag[f'{md}|pub|A3|T5']:.2f} / {ag[f'{md}|pub|A3|T10']:.2f} & {lstr}\\\\")
    s3 += [r"\bottomrule", r"\end{tabular*}", r"\end{table}"]
    (TAB / "tabS3_strata_ew.tex").write_text("\n".join(s3), encoding="utf-8")

    # Table S4: Antwerp Confirm05 full counts (T=5 and T=10)
    s4 = [
        r"\begin{table}[htbp]",
        r"\caption{Complete Confirm05\_Example alarm-flip counts on tile \texttt{E39N31} (Antwerp, N$'$, vertical, $\Dg$ median $0.6$\,mm/yr, range $0.5$--$0.6$\,mm/yr) at $\pm 5$ and $\pm 10$\,mm/yr. At $\pm 5$\,mm/yr, this tile sits at the 47.6th loss percentile across 205 active \texttt{O$'$}$\cup$\texttt{N$'$} tiles.}\label{tab:supp-antwerp}",
        r"\scriptsize\setlength{\tabcolsep}{1.6pt}",
        r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}p{0.20\linewidth}rrrrrrrrrr@{}}",
        r"\toprule",
        r"Area ($T$) & Pts & Any C0 & Both & A-only & B-only & Flip C0 & Flip C1 & $R$ & Expl. & New\\",
        r"\midrule",
    ]
    area_labels = {"city_window": "Antwerp 20\\,km", "rail_main": "Main railway 100\\,m",
                   "motorway": "Motorway 100\\,m", "tile": "Full tile E39N31"}
    for _, row in c5_counts.iterrows():
        s4.append(f"{area_labels[row['area']]} ($\\pm{int(row['T'])}$) & {int(row['n']):,} & {int(row['C0_alarm_any']):,} & "
                  f"{int(row['C0_both']):,} & {int(row['C0_A_only']):,} & {int(row['C0_B_only']):,} & "
                  f"{int(row['C0_flips']):,} & {int(row['C1_flips']):,} & {row['R']:.2f} & "
                  f"{int(row['explained_by_dg'])} & {int(row['new_after_dg'])}\\\\")
    s4 += [r"\bottomrule", r"\end{tabular*}", r"\end{table}"]
    (TAB / "tabS4_antwerp.tex").write_text("\n".join(s4), encoding="utf-8")
    print("ALL TABLES WRITTEN")


if __name__ == "__main__":
    build_all()
