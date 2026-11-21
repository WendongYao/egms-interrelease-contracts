"""Build all 6 main figures (Fig01-Fig06) and 5 supplementary figures (FigS1-FigS5)
for the 2026-11-21 revised JAG manuscript directly from local frozen contract results.

Nothing is recomputed or modified relative to registered contract outputs except:
- visual display jitter (+-0.035 mm/yr, seed=1) in Fig 2(a,b) matching draft_20261121;
- presentation-only redraw (z-order, marker sizes, figsize) of Confirm05_Example Antwerp E39N31
  maps (Fig06 & FigS5), with exact point counts asserted against example_counts.csv.
"""
import json
import os
import pathlib
import re
import numpy as np
import pandas as pd
import pyproj
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle
from shapely import STRtree, points
from shapely.geometry import LineString

ROOT = pathlib.Path(__file__).resolve().parents[2]
C01 = ROOT / "contracts" / "Confirm01_JAG"
C02 = ROOT / "contracts" / "Confirm02_Mechanism"
C03 = ROOT / "contracts" / "Confirm03_Proximity"
C04 = ROOT / "contracts" / "Confirm04_Robustness"
C05 = ROOT / "contracts" / "Confirm05_Example"
CACHE = ROOT / "manuscript" / "data_cache"
FIG = ROOT / "manuscript" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# Unified CVD-safe palette across ALL figures
C = {
    "O": "#9fc2ec",     # light blue for P1 pair O when contrasted against O'
    "O_map": "#2a78d6", # blue for shared 136 O/O' tiles in Fig 1
    "O'": "#2a78d6",    # blue for O' (P2)
    "N'": "#eb6834",    # vermilion/orange for N' (P2)
    "S''": "#1baf7a",   # green for S'' sparse supplement (P2)
}
INK, INK2, GRID, SHADE = "#0b0b0b", "#4f4e4b", "#d9d8d4", "#f2f1ee"
LAB = {"O": "O (P1)", "O'": "O′ (P2)", "N'": "N′ (P2)", "S''": "S″ (P2)"}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.0,
    "axes.titlesize": 8.8,
    "axes.labelsize": 8.0,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.edgecolor": INK2,
    "axes.labelcolor": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.45,
    "axes.axisbelow": True,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
    "figure.dpi": 150,
})


def save_fig(fig, name):
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight", dpi=250)
    plt.close(fig)
    print(f"wrote {name}.pdf/.png")


def tile_xy(tile_id):
    e, n = [int(v) for v in re.findall(r"\d+", tile_id)]
    return e * 1e5, n * 1e5


# ----------------------------------------------------------------------
# Fig 01: Study design, tile sets (with O central gap + S''), GNSS & flow
# ----------------------------------------------------------------------
def make_fig01():
    tr = pyproj.Transformer.from_crs(4326, 3035, always_xy=True)
    gj = json.load(open(CACHE / "ne_50m_admin_0_countries.geojson", encoding="utf-8"))

    mO = pd.read_csv(C01 / "data" / "manifest_O.csv")
    mN = pd.read_csv(C01 / "data" / "manifest_Nprime.csv")
    mS = pd.read_csv(C03 / "data" / "manifest_Spp.csv")

    bsta = pd.read_csv(C03 / "results" / "B_stations.csv")
    used_all = set(bsta["sta"])
    used_p2 = set(bsta[bsta["set"].isin(["O'", "N'"])]["sta"])
    st_meta = pd.concat([
        pd.read_csv(C01 / "data" / "stations_Oprime.csv")[["sta", "E", "N"]],
        pd.read_csv(C01 / "data" / "stations_Nprime.csv")[["sta", "E", "N"]],
    ]).drop_duplicates("sta")
    st_all = st_meta[st_meta["sta"].isin(used_all)].copy()
    assert len(used_p2) == 597 and len(st_all) == 618

    fig = plt.figure(figsize=(7.1, 4.7))
    # Map axis on the left
    ax = fig.add_axes([0.01, 0.03, 0.59, 0.94])
    ax.set_aspect("equal")
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

    for f in gj["features"]:
        g = f["geometry"]
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        for p in polys:
            ring = np.array(p[0])
            if ring[:, 1].max() < 25 or ring[:, 0].min() > 60 or ring[:, 0].max() < -35:
                continue
            x, y = tr.transform(ring[:, 0], ring[:, 1])
            ax.add_patch(Polygon(np.c_[x, y], closed=True, fc="#f1f0ec", ec="#b9b7b0", lw=0.35, zorder=1))

    # Highlight O's central-belt gap (27 <= k <= 35 -> N in [2.7e6, 3.6e6] m)
    ax.axhspan(2.7e6, 3.6e6, color="#e7e5de", alpha=0.45, zorder=2)
    ax.axhline(2.7e6, color=INK2, lw=0.6, ls="--", zorder=2)
    ax.axhline(3.6e6, color=INK2, lw=0.6, ls="--", zorder=2)
    ax.text(2.50e6, 3.56e6, "Central belt absent\nfrom O (27 ≤ k ≤ 35)", fontsize=7.0,
            color=INK2, va="top", ha="left", zorder=5,
            bbox=dict(fc="white", ec="none", alpha=0.72, pad=1.2))

    for t in mO["tile"]:
        x, y = tile_xy(t)
        ax.add_patch(Rectangle((x, y), 1e5, 1e5, fc=C["O_map"], ec="white", lw=0.3, alpha=0.88, zorder=3))
    for t, st in zip(mN["tile"], mN["stratum"]):
        x, y = tile_xy(t)
        ax.add_patch(Rectangle((x, y), 1e5, 1e5, fc=C["N'"],
                               ec=INK if st == "Al" else "white",
                               lw=0.7 if st == "Al" else 0.3, alpha=0.92, zorder=3))
    for t in mS["tile"]:
        x, y = tile_xy(t)
        ax.add_patch(Rectangle((x, y), 1e5, 1e5, fc=C["S''"], ec="white", lw=0.3, alpha=0.90, zorder=3))

    ax.scatter(st_all["E"], st_all["N"], s=2.6, c=INK, lw=0, zorder=4)
    ax.set_xlim(2.45e6, 6.55e6)
    ax.set_ylim(1.30e6, 5.45e6)
    ax.text(0.01, 0.99, "a  Tiles and GNSS stations (EPSG:3035)", transform=ax.transAxes,
            fontsize=8.8, fontweight="bold", va="top", color=INK)

    # Panel b: release timelines (top right)
    bx = fig.add_axes([0.65, 0.62, 0.34, 0.31])
    rels = [
        ("2019–2023", 2019, 2024, "XML v2.0\n(ETRF2000)", "#b5b3ac"),
        ("2020–2024", 2020, 2025, "XML v2.0\n(ETRF2000)", "#6e6c66"),
        ("2020–2025", 2020, 2026, "XML v3.0\n(ETRF2020)", INK),
    ]
    for i, (lab, a, b, gm, col) in enumerate(rels):
        y = 2 - i
        bx.barh(y, b - a, left=a, height=0.45, color=col, ec="none")
        bx.text(a + 0.08, y, lab, fontsize=7.2, color="white", va="center", fontweight="bold")
        bx.text(b + 0.12, y, gm, fontsize=7.0, color=INK2, va="center")
    bx.plot([2018.6, 2018.6], [1.0, 2.0], color=C["O'"], lw=1.4)
    bx.text(2018.45, 1.5, "P1", color=C["O'"], ha="right", va="center", fontsize=7.8, fontweight="bold")
    bx.plot([2018.1, 2018.1], [0.0, 1.0], color=C["N'"], lw=1.4)
    bx.text(2017.95, 0.5, "P2", color=C["N'"], ha="right", va="center", fontsize=7.8, fontweight="bold")
    bx.set_xlim(2017.3, 2029.6)
    bx.set_ylim(-0.6, 2.65)
    bx.set_yticks([])
    bx.set_xticks([2019, 2021, 2023, 2025])
    bx.tick_params(axis="x", labelsize=7.2)
    bx.spines["left"].set_visible(False)
    bx.grid(axis="y", visible=False)
    bx.set_xlabel("Release span (calendar years)", fontsize=7.6)
    bx.set_title("b  Release windows and pairs", loc="left", fontsize=8.8, fontweight="bold", color=INK)

    # Legend (middle right, below b)
    handles = [
        Rectangle((0, 0), 1, 1, fc=C["O_map"]),
        Rectangle((0, 0), 1, 1, fc=C["N'"]),
        Rectangle((0, 0), 1, 1, fc=C["N'"], ec=INK, lw=0.9),
        Rectangle((0, 0), 1, 1, fc=C["S''"]),
        Line2D([], [], marker="o", ls="", ms=3.2, color=INK),
    ]
    labels = [
        f"O / O′ tiles (n = {len(mO)}; both U and E on P2)",
        f"N′ stratified unseen tiles (n = {len(mN)}; U and E)",
        "N′ alpine stratum (n = 12, outlined)",
        f"S″ northern (Nt) suppl. (n = {len(mS)}; 15S/12M/8D, U)",
        f"NGL GNSS stations (618 total; 597 on P2)",
    ]
    fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.615, 0.555),
               handlelength=1.2, fontsize=7.2, labelspacing=0.45, frameon=False)

    # Panel c: compact pipeline & RSC strip (bottom right, replaces standalone box Fig 2)
    cx = fig.add_axes([0.605, 0.04, 0.39, 0.25])
    cx.set_xlim(0, 100)
    cx.set_ylim(0, 42)
    cx.axis("off")
    cx.text(0, 40, "c  Pair audit & Release Stability Check", fontsize=8.4, fontweight="bold", color=INK, va="top")

    steps = [
        (0.5, 20.5, "1. Match\n& XML", "#f6f5f2", INK2),
        (26.5, 20.5, "2. Frame\n−Δg (C1)", "#eef3fb", C["O'"]),
        (52.5, 20.5, "3. Plane\n(opt.)", "#f6f5f2", INK2),
        (78.5, 21.0, "4. Metrics\n& GHSL", "#fdf2ec", C["N'"]),
    ]
    for sx, sw, txt, fc, ec in steps:
        cx.add_patch(FancyBboxPatch((sx, 13.5), sw, 18.0, boxstyle="round,pad=0.25,rounding_size=1.4",
                                    fc=fc, ec=ec, lw=0.75))
        cx.text(sx + sw / 2.0, 22.5, txt, ha="center", va="center", fontsize=6.9, color=INK, linespacing=1.22)
    for ax_x0, ax_x1 in ((21.4, 26.1), (47.4, 52.1), (73.4, 78.1)):
        cx.annotate("", xy=(ax_x1, 22.5), xytext=(ax_x0, 22.5),
                    arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INK2, mutation_scale=7))
    cx.add_patch(FancyBboxPatch((0.5, 0.8), 99.0, 9.2, boxstyle="round,pad=0.2,rounding_size=1.0",
                                fc="#f1f0ec", ec="none"))
    cx.text(50.0, 5.4, "Contracts (Confirm01–05) logged pre-run; all gates reported",
            ha="center", va="center", fontsize=6.6, color=INK2)

    save_fig(fig, "Fig01_study_design")


# ----------------------------------------------------------------------
# Fig 02: Reference-frame term decomposition & removal (U, E, S0 -> S1)
# ----------------------------------------------------------------------
def make_fig02():
    bu_o = pd.read_csv(C02 / "results" / "mech" / "by_tile_O'.csv")
    bu_n = pd.read_csv(C02 / "results" / "mech" / "by_tile_N'.csv")
    bu_o = bu_o[bu_o.included == True]
    bu_n = bu_n[bu_n.included == True]
    du = pd.read_csv(C03 / "results" / "D_tile_U.csv")
    de = pd.read_csv(C03 / "results" / "D_tile_E.csv")
    r1 = pd.read_csv(C04 / "results" / "R1.csv").set_index("set")

    rng = np.random.default_rng(1)
    jit = lambda n: rng.uniform(-0.035, 0.035, n)

    fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.75),
                            gridspec_kw=dict(width_ratios=[1, 1, 1.15], wspace=0.33))

    panels = [
        (axs[0], "a  Vertical (U)", (-1.0, 2.0),
         [("O'", bu_o), ("N'", bu_n), ("S''", du[du.set == "S''"])]),
        (axs[1], "b  East–west (E)", (-2.0, 3.0),
         [("O'", de[de.set == "O'E"]), ("N'", de[de.set == "N'E"])]),
    ]
    for ax, title, lim, series in panels:
        for name, df in series:
            ax.scatter(df["med_dg"] + jit(len(df)), df["med_dv_pub"] + jit(len(df)),
                       s=9, color=C[name], alpha=0.78, lw=0, label=LAB[name])
        ax.plot(lim, lim, color=INK2, lw=0.8, ls="--")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel("Tile median Δg (mm/yr)")
        ax.set_title(title, loc="left", fontsize=8.8, fontweight="bold", color=INK)
    axs[0].set_ylabel("Tile median Δv, published (mm/yr)")
    axs[0].legend(fontsize=7.2, loc="upper left", handletextpad=0.3)

    ax = axs[2]
    keys = [("O'U", "O′ U"), ("N'U", "N′ U"), ("O'E", "O′ E"), ("N'E", "N′ E"), ("S''", "S″ U")]
    yy = np.arange(len(keys))[::-1]
    for y, (k, lbl) in zip(yy, keys):
        row = r1.loc[k]
        ax.plot([row["S1_hd"], row["S0_hd"]], [y, y], color=GRID, lw=2.2)
        ax.scatter(row["S0_hd"], y, s=26, color=INK2, zorder=3)
        ax.scatter(row["S1_hd"], y, s=26, color=C["N'"], zorder=3)
        ax.text(row["S0_hd"] + 0.03, y + 0.14, f"{row['S0_hd']:.2f}", fontsize=7.0, color=INK2)
        s1_txt = "0.00" if abs(row["S1_hd"]) < 0.005 else f"{row['S1_hd']:.2f}"
        ax.text(max(row["S1_hd"], 0.0) + 0.035, y - 0.34, s1_txt, fontsize=7.0, color=C["N'"])
    ax.set_yticks(yy)
    ax.set_yticklabels([t for _, t in keys])
    ax.axvline(0, color=INK2, lw=0.65)
    ax.set_xlim(-0.14, 0.88)
    ax.set_ylim(-0.6, len(keys) - 0.35)
    ax.scatter([], [], s=24, color=INK2, label="before (S0)")
    ax.scatter([], [], s=24, color=C["N'"], label="after −Δg (S1)")
    ax.legend(fontsize=7.2, loc="upper center", bbox_to_anchor=(0.48, -0.25), ncol=2, handletextpad=0.3)
    ax.set_xlabel("Systematic offset (Harrell–Davis, mm/yr)")
    ax.set_title("c  Offset before / after", loc="left", fontsize=8.8, fontweight="bold", color=INK)

    fig.tight_layout()
    save_fig(fig, "Fig02_frame_term")


# ----------------------------------------------------------------------
# Fig 03: Decision instability (C0 -> C2) & frame share s_g
# ----------------------------------------------------------------------
def make_fig03():
    s2all = json.load(open(C02 / "results" / "mech" / "summary.json", encoding="utf-8"))
    s2_shares = s2all["M5"]["shares"]
    ag2 = s2all["aggregates"]
    sd = json.load(open(C03 / "results" / "summary_D.json", encoding="utf-8"))["D1"]

    rows = [
        ("A1", "T0", "Point flips, sign (v < 0)"),
        ("A1", "T05", "Point flips, ±0.5 mm/yr"),
        ("A1", "T5", "Point flips, ±5 mm/yr"),
        ("A1", "T10", "Point flips, ±10 mm/yr"),
        ("A2", "T5", "Hotspot loss, ±5 mm/yr"),
        ("A2", "T10", "Hotspot loss, ±10 mm/yr"),
        ("A3", "T5", "Built-up disagreement, ±5"),
        ("A3", "T10", "Built-up disagreement, ±10"),
    ]
    series = [
        ("O'", "U", C["O'"], "o", "O′ U"),
        ("N'", "U", C["N'"], "o", "N′ U"),
        ("O'", "E", C["O'"], "^", "O′ E"),
        ("N'", "E", C["N'"], "^", "N′ E"),
    ]

    fig, (ax0, ax) = plt.subplots(1, 2, figsize=(7.1, 3.55), sharey=True,
                                  gridspec_kw=dict(width_ratios=[1.0, 1.12], wspace=0.10))
    yy = np.arange(len(rows))[::-1]
    for j, (n, comp, col, mk, lbl) in enumerate(series):
        off = (j - 1.5) * 0.17
        for y, (met, th, _) in zip(yy, rows):
            if comp == "E" and th == "T0":
                continue
            if comp == "U":
                d0 = ag2[f"{n}|C0|pub|{met}|{th}"]
                d2 = ag2[f"{n}|C2|pub|{met}|{th}"]
                s = s2_shares.get(f"{n}|pub|{met}|{th}")
            else:
                ag = sd[f"{n}E"]["metrics"]["aggregates"]
                d0 = ag[f"C0|pub|{met}|{th}"]
                d2 = ag[f"C2|pub|{met}|{th}"]
                s = sd[f"{n}E"]["metrics"]["shares"].get(f"pub|{met}|{th}")
            fc = col if comp == "U" else "white"
            ax0.plot([d0, d2], [y + off, y + off], color=col, lw=0.9, alpha=0.65)
            ax0.plot(d0, y + off, mk, ms=4.0, color=col, mfc="white", mew=0.9)
            ax0.plot(d2, y + off, mk, ms=4.0, color=col, mfc=col, mew=0.9)
            if s and "s_g" in s:
                lo_e = max(0.0, s["s_g"] - s["s_g_ci"][0])
                hi_e = max(0.0, s["s_g_ci"][1] - s["s_g"])
                ax.errorbar(s["s_g"], y + off, xerr=[[lo_e], [hi_e]],
                            fmt=mk, ms=4.0, color=col, mfc=fc, mew=0.9, elinewidth=0.9, capsize=0)

    ax0.set_xscale("log")
    ax0.set_xlim(3e-4, 1.2)
    ax0.set_xlabel("Disagreement (log scale)\nopen: as published · filled: −Δg & tile plane")
    ax0.set_yticks(yy)
    ax0.set_yticklabels([r[2] for r in rows])
    ax0.set_title("a  Disagreement before / after", loc="left", fontsize=8.8, fontweight="bold", color=INK)

    ax.axvline(0, color=INK2, lw=0.7)
    ax.set_xlabel(r"Frame share $s_g$" + "\n(95% CI, tile bootstrap)")
    ax.set_xlim(-0.65, 0.85)
    for n, comp, col, mk, lbl in series:
        ax.plot([], [], mk, ms=4.5, color=col, mfc=col if comp == "U" else "white", mew=0.9, label=lbl)
    ax.legend(fontsize=7.3, loc="lower right", ncol=2, columnspacing=0.9, handletextpad=0.3)
    ax.set_title(r"b  Frame contribution share $s_g$", loc="left", fontsize=8.8, fontweight="bold", color=INK)

    fig.tight_layout()
    save_fig(fig, "Fig03_decisions")


# ----------------------------------------------------------------------
# Fig 04: Threshold sensitivity & boundary proximity
# ----------------------------------------------------------------------
def make_fig04():
    cur = pd.read_csv(C04 / "results" / "R3_curves.csv")
    ab = pd.read_csv(C03 / "results" / "A_bins.csv")
    J = {"O": ("O", "C0"), "O'": ("O'U", "C1"), "N'": ("N'U", "C1")}

    fig, axs = plt.subplots(1, 4, figsize=(7.1, 2.50), gridspec_kw=dict(wspace=0.32))
    for n, (k, md) in J.items():
        d = cur[cur["set"] == k].sort_values("T")
        col = C[n]
        axs[0].plot(d["T"], d["loss_pooled"], color=col, lw=1.7, label=LAB[n])
        axs[1].plot(d["T"], d["rho1_loss"], color=col, lw=1.7)
        axs[1].plot(d["T"], d["rho1_gain"], color=col, lw=1.1, ls="--")
        axs[2].plot(d["T"], d["rho2"], color=col, lw=1.7)
        b = ab[(ab["set"] == n) & (ab["var"] == "pub") & (ab["mode"] == md) & (ab["T"] == 5) & (ab["dir"] == "loss")]
        axs[3].plot(range(len(b)), b["rate"].values, "o-", ms=3.5, color=col, lw=1.3)

    for a in axs[:3]:
        a.axvspan(3, 12, color=SHADE, zorder=0)
        a.set_xlabel("Threshold T (mm/yr)")
        a.set_xticks([2, 5, 10, 15])
    axs[1].axhline(1, color=INK2, lw=0.65)
    axs[2].axhline(1, color=INK2, lw=0.65)
    axs[0].set_ylabel("Fraction")
    axs[0].set_title("a  Pooled loss rate", loc="left", fontsize=8.4, fontweight="bold", color=INK)
    axs[1].set_title("b  ρ₁ vs random sign", loc="left", fontsize=8.4, fontweight="bold", color=INK)
    axs[2].set_title("c  ρ₂ vs typical-point", loc="left", fontsize=8.4, fontweight="bold", color=INK)
    axs[3].set_title("d  Loss by margin (±5)", loc="left", fontsize=8.4, fontweight="bold", color=INK)

    axs[0].legend(fontsize=7.2, loc="lower right", handlelength=1.4)
    axs[1].plot([], [], color=INK2, lw=1.7, label="loss")
    axs[1].plot([], [], color=INK2, lw=1.1, ls="--", label="gain")
    axs[1].legend(fontsize=7.2, loc="lower right", handlelength=1.4)
    axs[3].set_xticks(range(5))
    axs[3].set_xticklabels(["0–0.5", "0.5–1", "1–2", "2–5", "≥5"], rotation=32, ha="right", fontsize=7.2)
    axs[3].set_xlabel(r"$|v_A| - 5$ (mm/yr)")

    fig.tight_layout()
    save_fig(fig, "Fig04_threshold")


# ----------------------------------------------------------------------
# Fig 05: GNSS offset decomposition & station-level R_dec (300 km block CI)
# ----------------------------------------------------------------------
def make_fig05():
    b2 = json.load(open(C03 / "results" / "summary.json", encoding="utf-8"))["B2"]
    r2 = pd.read_csv(C04 / "results" / "R2.csv").set_index("key")

    fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.75),
                            gridspec_kw=dict(width_ratios=[1.38, 1.0], wspace=0.28))

    ax = axs[0]
    items = [
        ("P1 (O)", "2019–23", b2["P1"]["c_2019_23"], "#9e9c96"),
        ("P1 (O)", "2020–24", b2["P1"]["c_2020_24_P1"], INK2),
    ]
    for n in ("O'", "N'"):
        tag = f"P2 {n}"
        items += [
            (tag, "2020–24", b2[n]["c_2020_24_P2"], INK2),
            (tag, "2020–25", b2[n]["c_2020_25"], INK),
            (tag, "2020–25 −Δg", b2[n]["c_2020_25_C1"], C["N'"]),
        ]
    xs, gap = [], 0.0
    for i, (_, _, v, col) in enumerate(items):
        if i in (2, 5):
            gap += 0.65
        xi = i + gap
        xs.append(xi)
        ax.errorbar(xi, v["est"], yerr=[[v["est"] - v["ci"][0]], [v["ci"][1] - v["est"]]],
                    fmt="o", ms=4.5, color=col, elinewidth=1.1, capsize=0)
    ax.set_xticks(xs)
    ax.set_xticklabels([r for _, r, _, _ in items], rotation=32, ha="right", fontsize=7.2)
    for xi, t in ((0.5, "P1 (O)"), (3.65, "P2 O′"), (7.30, "P2 N′")):
        ax.text(xi, 0.06, t, ha="center", fontsize=7.6, color=INK, fontweight="bold")
    ax.axhline(0, color=INK2, lw=0.65)
    ax.set_ylabel("EGMS − GNSS offset c (mm/yr)")
    ax.set_ylim(-1.45, 0.22)
    ax.set_title("a  Per-release offset at GNSS stations", loc="left", fontsize=8.8, fontweight="bold", color=INK)

    ax = axs[1]
    ks = [("O'", "T0", "O′ sign"), ("O'", "T05", "O′ ±0.5"),
          ("N'", "T0", "N′ sign"), ("N'", "T05", "N′ ±0.5")]
    yy = np.arange(len(ks))[::-1]
    for y, (n, th, _) in zip(yy, ks):
        for md, col, dy in (("C0", "#9e9c96", 0.15), ("C1", C["N'"], -0.15)):
            r = r2.loc[f"Rdec|{n}|{md}|{th}"]
            ax.errorbar(r["est"], y + dy, xerr=[[r["est"] - r["lo300"]], [r["hi300"] - r["est"]]],
                        fmt="o", ms=4.5, color=col, elinewidth=1.1, capsize=0)
    ax.axvline(1, color=INK2, lw=0.7)
    ax.set_yticks(yy)
    ax.set_yticklabels([t for _, _, t in ks])
    ax.set_ylim(-0.55, 4.15)
    ax.plot([], [], "o", ms=4.5, color="#9e9c96", label="as pub. (C0)")
    ax.plot([], [], "o", ms=4.5, color=C["N'"], label="after −Δg (C1)")
    ax.legend(fontsize=7.0, loc="upper right", ncol=2, columnspacing=0.9, handletextpad=0.3)
    ax.set_xlabel(r"Station flip ratio $R_{\mathrm{dec}}$ (300 km block CI)")
    ax.set_title("b  Release flips vs EGMS–GNSS flips", loc="left", fontsize=8.8, fontweight="bold", color=INK)

    fig.tight_layout()
    save_fig(fig, "Fig05_gnss")


# ----------------------------------------------------------------------
# Fig 06 & FigS5: Confirm05 Antwerp E39N31 example (redrawn, OSM over grey dots)
# ----------------------------------------------------------------------
def _load_e39n31():
    z = np.load(C01 / "data" / "fulltile_Nprime" / "E39N31.npz")
    g = np.load(C02 / "data" / "g_Nprime" / "E39N31.npz")
    sel = json.load(open(C05 / "results" / "selection.json", encoding="utf-8"))
    summ = json.load(open(C05 / "results" / "summary.json", encoding="utf-8"))
    exp = pd.read_csv(C05 / "results" / "example_counts.csv")
    assert len(z["E"]) == len(g["gA"]) == 353140
    d = dict(E=z["E"].astype(float), N=z["N"].astype(float),
             vA=z["vA_pub"].astype(float), vB=z["vB_pub"].astype(float),
             gA=g["gA"].astype(float), gB=g["gB"].astype(float))

    js = json.load(open(C05 / "data" / "osm" / "E39N31_overpass.json", encoding="utf-8"))
    tr = pyproj.Transformer.from_crs(4326, 3035, always_xy=True)
    rail, mway = [], []
    for el in js["elements"]:
        tg = el.get("tags", {})
        if el["type"] == "way" and "geometry" in el and len(el["geometry"]) >= 2:
            x, y = tr.transform([p["lon"] for p in el["geometry"]], [p["lat"] for p in el["geometry"]])
            ln = LineString(list(zip(x, y)))
            if tg.get("railway") == "rail" and tg.get("usage") == "main":
                rail.append(ln)
            elif tg.get("highway") == "motorway":
                mway.append(ln)
    return d, rail, mway, sel, summ, exp


def _near_lines(E, N, lines, buf=100.0):
    tree = STRtree(lines)
    pts = points(np.c_[E, N])
    idx = tree.query(pts, predicate="dwithin", distance=buf)
    m = np.zeros(len(E), bool)
    m[np.unique(idx[0])] = True
    return m


def _states(a, b, T):
    A = np.abs(a) >= T
    B = np.abs(b) >= T
    return np.select([A & B, A & ~B, ~A & B], [0, 1, 2], 3)


def make_antwerp_fig(T, out_name, d, rail, mway, sel, exp, masks):
    col = {0: C["O'"], 1: C["N'"], 2: C["S''"], 3: "#cccccc"}
    lab = {0: "alarm in both", 1: "alarm only in 2020–2024", 2: "alarm only in 2020–2025", 3: "no alarm"}
    dg = d["gB"] - d["gA"]
    ce, cn = sel["window_centre_E"], sel["window_centre_N"]
    h = sel["window_side_m"] / 2.0
    e0, n0 = sel["tile_origin_E"], sel["tile_origin_N"]

    titles = {
        "city_window": "Antwerp (20 km window)",
        "rail_main": "Main railways (100 m buffer)",
        "motorway": "Motorways (100 m buffer)",
    }

    fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.95), constrained_layout=True)
    for j, area in enumerate(("city_window", "rail_main", "motorway")):
        m = masks[area]
        E_km = d["E"][m] / 1000.0
        N_km = d["N"][m] / 1000.0
        s0 = _states(d["vA"][m], d["vB"][m], T)
        s1 = _states(d["vA"][m], d["vB"][m] - dg[m], T)

        # Verify exact match to frozen Confirm05 CSV counts
        ref_row = exp[(exp["area"] == area) & (exp["T"] == float(T))].iloc[0]
        assert int(m.sum()) == int(ref_row["n"]), f"{area} T{T} n mismatch"
        assert int(np.isin(s0, (1, 2)).sum()) == int(ref_row["C0_flips"]), f"{area} T{T} C0 flips"
        assert int(np.isin(s1, (1, 2)).sum()) == int(ref_row["C1_flips"]), f"{area} T{T} C1 flips"

        lines_to_plot = (rail if area != "motorway" else []) + (mway if area != "rail_main" else [])

        for i, (s, tag) in enumerate(((s0, "C0: as published"), (s1, "C1: 2020–2025 −Δg"))):
            ax = axs[i, j]
            ax.grid(False)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(True)
                ax.spines[sp].set_color(INK2)
                ax.spines[sp].set_linewidth(0.5)

            # 1) Background no-alarm points BELOW OSM lines
            q_bg = (s == 3)
            ax.scatter(E_km[q_bg], N_km[q_bg], s=0.35 if area == "city_window" else 0.45,
                       c=col[3], alpha=0.45, lw=0, zorder=1, rasterized=True)

            # 2) OSM linear infrastructure ABOVE no-alarm points (LineCollection for speed)
            segs = [np.column_stack(ln.xy) / 1000.0 for ln in lines_to_plot]
            if segs:
                lc = LineCollection(segs, colors="#3b3a38", linewidths=0.65, alpha=0.88, zorder=2)
                ax.add_collection(lc)

            # 3) Persistent and flipped alarm points ON TOP
            for k, sz, zo in ((0, 3.2, 3), (1, 5.5, 4), (2, 5.5, 4)):
                q = (s == k)
                if q.any():
                    ax.scatter(E_km[q], N_km[q], s=sz, c=col[k], lw=0, zorder=zo)

            if area == "city_window":
                ax.set_xlim((ce - h) / 1000.0, (ce + h) / 1000.0)
                ax.set_ylim((cn - h) / 1000.0, (cn + h) / 1000.0)
            else:
                ax.set_xlim(e0 / 1000.0, e0 / 1000.0 + 100.0)
                ax.set_ylim(n0 / 1000.0, n0 / 1000.0 + 100.0)
            ax.set_aspect("equal")
            ax.tick_params(labelsize=7.0)
            nf = int(np.isin(s, (1, 2)).sum())
            ax.text(0.03, 0.96, f"flips: {nf:,}", transform=ax.transAxes, va="top", fontsize=7.5,
                    fontweight="bold", color=INK,
                    bbox=dict(facecolor="white", edgecolor=GRID, lw=0.4, alpha=0.88, pad=1.8), zorder=5)
            ax.set_title(f"{titles[area]}\n{tag}", fontsize=8.0, color=INK)
            if j == 0:
                ax.set_ylabel("Northing (km)", fontsize=7.5)
            if i == 1:
                ax.set_xlabel("Easting (km, EPSG:3035)", fontsize=7.5)

    handles = [Line2D([], [], marker="o", ls="", color=col[k], label=lab[k], markersize=4.8) for k in (0, 1, 2, 3)]
    handles.append(Line2D([], [], color="#3b3a38", lw=1.0, label="OSM rail / motorway"))
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=7.2,
               frameon=False, bbox_to_anchor=(0.5, -0.055), handletextpad=0.35, columnspacing=1.0)
    save_fig(fig, out_name)


# ----------------------------------------------------------------------
# Supplementary figures (FigS1-FigS5)
# ----------------------------------------------------------------------
def make_figs_supp(d, rail, mway, sel, exp, masks):
    # Fig S1: Detailed study workflow & RSC diagram
    fig, ax = plt.subplots(figsize=(7.1, 4.0))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 56)
    ax.axis("off")

    def box(x, y, w, h, title, lines, fc="#ffffff", ec=INK2):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.2,rounding_size=1.2",
                                    fc=fc, ec=ec, lw=0.75))
        ax.text(x + 1.0, y + h - 1.1, title, fontsize=7.5, fontweight="bold", va="top", color=INK)
        ax.text(x + 1.0, y + h - 4.3, "\n".join(lines), fontsize=7.0, va="top", color=INK2, linespacing=1.35)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INK2, mutation_scale=7.5))

    box(0.5, 29, 20.5, 25, "Inputs", ["EGMS Ortho U & E", "(P1, P2; 243 tiles)", "XML GNSS model", "NGL GNSS (618 stns)", "GHSL R2023A grids"])
    box(23, 29, 20, 25, "Alignment", ["Match 100 m points", "As-published and", "common-window OLS", "Per-point Δg from", "gnss_velocity_{u,e}"])
    box(45.5, 42, 24.5, 12, "Point differences", ["Δv vs formal errors,", "density, GNSS (F1–F11)"], fc="#f6f5f2")
    box(45.5, 29, 24.5, 12, "Screening decisions", ["0, ±0.5, ±5, ±10 mm/yr", "hotspots, GHSL (A1–A4)"], fc="#f6f5f2")
    box(72.5, 29, 27, 25, "Attribution", ["C0: as published", "C1: C0 − Δg (frame)", "C2: C1 − tile plane", "Shares s_g, s_ramp, s_local", "Proximity & blocks"])
    box(0.5, 7.5, 69.5, 17.5, "Release Stability Check (RSC)",
        ["1 Metadata (XML version) → 2 Frame (−Δg / constant) →",
         "3 Optional tile plane (C2 / C2′) → 4 Report (A1–A3, hotspots)",
         "Checks U1–U3: bit-identical (162 cols), injection, Δg≡0 control"],
        fc="#eef3fb", ec=C["O'"])
    box(72.5, 7.5, 27, 17.5, "Decision rules",
        ["Near 0: remove Δg first", "±5: 37–62% tile-median", "(16–22% pooled) C0 loss"])
    arrow(20.7, 41.5, 22.8, 41.5)
    arrow(43.2, 48, 45.3, 48)
    arrow(43.2, 35, 45.3, 35)
    arrow(70.2, 41.5, 72.3, 41.5)
    arrow(86, 28.8, 86, 25.2)
    arrow(70.2, 16.2, 72.3, 16.2)
    arrow(33, 28.8, 33, 25.2)
    ax.add_patch(FancyBboxPatch((0.5, 0.8), 99, 5.0, boxstyle="round,pad=0.1,rounding_size=0.8",
                                fc="#f1f0ec", ec="none"))
    ax.text(50, 3.3, "SHA-256 hashes logged locally pre-run; public archiving post-dates analysis (not third-party timestamped).",
            ha="center", va="center", fontsize=6.7, color=INK)
    save_fig(fig, "FigS_workflow")

    # Fig S2: Uncertainty diagnostics (F4, F5, F6) with unified palette
    s_o = json.load(open(C01 / "results" / "part1_Oprime__summary.json", encoding="utf-8"))
    s_n = json.load(open(C01 / "results" / "part1_Nprime__summary.json", encoding="utf-8"))
    b_o = pd.read_csv(C01 / "results" / "part1_Oprime__by_point.csv", usecols=["kind", "mean_velocity_std"])
    b_n = pd.read_csv(C01 / "results" / "part1_Nprime__by_point.csv", usecols=["kind", "mean_velocity_std"])
    S_dict = {"O'": (s_o, b_o), "N'": (s_n, b_n)}

    fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.7), gridspec_kw=dict(wspace=0.38, width_ratios=[1, 1.15, 1]))
    cl = ["dense", "mid", "sparse"]
    axs[0].axhspan(0.92, 0.97, color=SHADE, lw=0)
    axs[0].text(2.42, 0.945, "target\n0.92–0.97", fontsize=7.0, color=INK2, va="center", ha="right")
    for j, (n, (s, _)) in enumerate(S_dict.items()):
        v = [s["F4"]["coverage"]["by_class"][c]["coverage"] for c in cl]
        axs[0].plot(np.arange(3) + (j - 0.5) * 0.18, v, "o", color=C[n], ms=5.0, label=LAB[n])
    axs[0].set_xticks(range(3))
    axs[0].set_xticklabels(cl)
    axs[0].set_ylim(0.78, 1.0)
    axs[0].set_xlim(-0.5, 2.5)
    axs[0].set_ylabel("Coverage of inflated 2σ interval")
    axs[0].axhline(0.88, color=INK2, lw=0.6, ls="--")
    axs[0].text(-0.45, 0.883, "floor 0.88", fontsize=7.0, color=INK2)
    axs[0].legend(loc="lower right", fontsize=7.2)
    axs[0].grid(axis="x", visible=False)
    axs[0].set_title("a  Coverage (P1 B0)", loc="left", fontsize=8.5, fontweight="bold", color=INK)

    for j, (n, (s, b)) in enumerate(S_dict.items()):
        v = b[b["kind"] == "rand"]["mean_velocity_std"].round(1)
        vc = v.value_counts(normalize=True).sort_index()
        axs[1].bar(vc.index + (j - 0.5) * 0.035, vc.values, width=0.033, color=C[n],
                   label=f"{LAB[n]}: {s['F5']['n_distinct_mvs']} values")
    axs[1].set_xlim(-0.07, 0.47)
    axs[1].set_ylim(0, 0.68)
    axs[1].set_xticks([0, 0.1, 0.2, 0.3, 0.4])
    axs[1].set_xlabel("mean_velocity_std (mm/yr)")
    axs[1].set_ylabel("Share of points")
    axs[1].legend(fontsize=7.2, loc="upper right")
    axs[1].grid(axis="x", visible=False)
    axs[1].set_title("b  Product std values", loc="left", fontsize=8.5, fontweight="bold", color=INK)

    for j, (n, (s, _)) in enumerate(S_dict.items()):
        v = [s["F6"]["by_class"][c]["r_median"] for c in cl]
        axs[2].plot(np.arange(3) + (j - 0.5) * 0.18, v, "o", color=C[n], ms=5.0)
        axs[2].plot([2.75 + (j - 0.5) * 0.18], [s["F6"]["r_median"]], "D", color=C[n], ms=4.5)
        axs[2].plot([2.75 + (j - 0.5) * 0.18] * 2, [s["F6"]["r_ci_lo"], s["F6"]["r_ci_hi"]], color=C[n], lw=1.0)
    axs[2].axhline(1, color=INK2, lw=0.65)
    axs[2].axhline(0.8, color=INK2, lw=0.6, ls="--")
    axs[2].set_xticks([0, 1, 2, 2.75])
    axs[2].set_xticklabels(cl + ["all"], fontsize=7.2)
    axs[2].set_ylim(0, 1.1)
    axs[2].set_ylabel("Short / long contrast ratio")
    axs[2].grid(axis="x", visible=False)
    axs[2].set_title("c  Spatial contrast (F6)", loc="left", fontsize=8.5, fontweight="bold", color=INK)
    fig.tight_layout()
    save_fig(fig, "FigS_uncert")

    # Fig S3: E32N16 hotspot example (rename copy of existing Fig04_example_tile)
    for ext in ("pdf", "png"):
        src = FIG / f"Fig04_example_tile.{ext}"
        dst = FIG / f"FigS_E32N16.{ext}"
        if src.exists():
            dst.write_bytes(src.read_bytes())
    print("copied FigS_E32N16.pdf/.png from Fig04_example_tile")

    # Fig S4: Density step on N' U S'' (C0 vs C1, common window)
    du = pd.read_csv(C03 / "results" / "D_tile_U.csv")
    U = du[du["included"] == True]
    fig, axs = plt.subplots(1, 2, figsize=(6.8, 2.55), sharey=True)
    for ax, md, title in ((axs[0], "C0", "a  As published (C0)"), (axs[1], "C1", "b  Frame-harmonized (C1, −Δg)")):
        for n in ("N'", "S''"):
            sub = U[U["set"] == n]
            ax.scatter(sub["n"], sub[f"{md}__mstar_cw"], s=11, color=C[n], alpha=0.82, lw=0, label=LAB[n])
        ax.axvspan(0, 130, color=SHADE, zorder=0)
        ax.axvline(200, color=GRID, lw=0.9, ls="--")
        ax.set_ylim(0, 3.05)
        ax.set_xlabel("Common-window epochs, n")
        ax.set_title(title, loc="left", fontsize=8.6, fontweight="bold", color=INK)
    axs[0].set_ylabel("Tile median |Δv| (mm/yr)")
    axs[0].legend(fontsize=7.3, loc="upper right")
    fig.tight_layout()
    save_fig(fig, "FigS_density_step")

    # Fig S5: Antwerp T10 companion map
    make_antwerp_fig(10.0, "FigS_antwerp_T10", d, rail, mway, sel, exp, masks)


def main():
    make_fig01()
    make_fig02()
    make_fig03()
    make_fig04()
    make_fig05()

    d, rail, mway, sel, summ, exp = _load_e39n31()
    ce, cn = sel["window_centre_E"], sel["window_centre_N"]
    h = sel["window_side_m"] / 2.0
    masks = {
        "city_window": (np.abs(d["E"] - ce) <= h) & (np.abs(d["N"] - cn) <= h),
        "rail_main": _near_lines(d["E"], d["N"], rail, 100.0),
        "motorway": _near_lines(d["E"], d["N"], mway, 100.0),
    }
    make_antwerp_fig(5.0, "Fig06_antwerp_T5", d, rail, mway, sel, exp, masks)
    make_figs_supp(d, rail, mway, sel, exp, masks)
    print("ALL FIGURES COMPLETE")


if __name__ == "__main__":
    main()
