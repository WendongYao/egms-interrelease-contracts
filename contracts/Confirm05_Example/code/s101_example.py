"""Confirm05_Example: one city window and two kinds of linear infrastructure (OSM) in one tile.

Stages (contract CONTRACT_CN.md, SHA-256 fa8f0f9c...0446):
  select   -> results/selection.json   (tile with max GHSL 2020 population at matched points; window centre)
  analyze  -> results/example_counts.csv, results/summary.json, results/Fig_example_T5.png, results/Fig_example_T10.png
  smoke    -> loads 2 tiles and an OSM file if present; prints field names and shapes only (no readings)

Inputs (cloud, unchanged outputs of earlier contracts):
  Confirm01 full-tile matches  /home/claude/c01/run/fulltile_{Oprime,Nprime}/<tile>.npz  (E, N, vA_pub, vB_pub)
  Confirm02 per-point GNSS U   /home/claude/c02/data/g_{Oprime,Nprime}/<tile>.npz       (gA, gB)
  GHSL 2020 POP (100 m)        /home/claude/c01/ghsl, read with Confirm01 s89 exposure() (copied verbatim in
                               Confirm02 code/rsc/p2lib_exposure.py)
  Confirm03 A_tile.csv         tile-level T5 C0 pooled loss for the percentile (contract section 4)
  OSM Overpass JSON            data/osm/<tile>_overpass.json (fetched on the VM by s100_osm_fetch.py)
"""
import os, sys, json, glob, hashlib
import numpy as np

C01 = "/home/claude/c01/run"; C02 = "/home/claude/c02/data"; C03 = "/home/claude/c03/results"
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results"); OSM = os.path.join(ROOT, "data", "osm")
os.environ.setdefault("GHSL_DIR", "/home/claude/c01/ghsl")
sys.path.insert(0, "/home/claude/c02/code/rsc")
SETS = {"Oprime": "O'", "Nprime": "N'"}
WIN = 20000.0       # city window side, m (contract 1)
BUF = 100.0         # infrastructure buffer, m (contract 1)
TS = (5.0, 10.0)    # T5 main, T10 report only (contract 3)


def tiles():
    out = []
    for s in SETS:
        for f in sorted(glob.glob(f"{C01}/fulltile_{s}/*.npz")):
            out.append((s, os.path.basename(f)[:-4]))
    return out


def load(s, t):
    z = np.load(f"{C01}/fulltile_{s}/{t}.npz"); g = np.load(f"{C02}/g_{s}/{t}.npz")
    assert len(g["gA"]) == len(z["E"])
    return dict(E=z["E"].astype(float), N=z["N"].astype(float), vA=z["vA_pub"].astype(float),
                vB=z["vB_pub"].astype(float), gA=g["gA"].astype(float), gB=g["gB"].astype(float))


def wmedian(x, w):
    o = np.argsort(x); x, w = x[o], w[o]; c = np.cumsum(w)
    return float(x[np.searchsorted(c, 0.5 * c[-1])])


def tile_box_ll(t):
    import pyproj
    e0 = int(t[1:t.index("N")]) * 100000; n0 = int(t[t.index("N") + 1:]) * 100000
    tr = pyproj.Transformer.from_crs(3035, 4326, always_xy=True)
    xs = np.r_[np.linspace(e0, e0 + 1e5, 21), np.full(21, e0 + 1e5), np.linspace(e0 + 1e5, e0, 21), np.full(21, e0)]
    ys = np.r_[np.full(21, n0), np.linspace(n0, n0 + 1e5, 21), np.full(21, n0 + 1e5), np.linspace(n0 + 1e5, n0, 21)]
    lo, la = tr.transform(xs, ys)
    return dict(south=float(la.min()), west=float(lo.min()), north=float(la.max()), east=float(lo.max())), (e0, n0)


def stage_select():
    import p2lib_exposure as px
    rows = []
    for s, t in tiles():
        d = load(s, t); ex = px.exposure(d["E"], d["N"])
        # raster coverage: points whose GHSL zip is absent get 0; count them so a missing raster is visible
        x, y = px.TR.transform(d["E"], d["N"])
        r = ((9000000 - y) // 1e6).astype(int) + 1; c = ((x + 18041000) // 1e6).astype(int) + 1
        miss = sum(int(((r == rr) & (c == cc)).sum()) for rr, cc in set(zip(r.tolist(), c.tolist()))
                   if px.ghsl_path("pop", rr, cc) is None)
        rows.append(dict(set=SETS[s], key=s, tile=t, n=len(d["E"]), pop_sum=float(ex["pop"].sum()), n_no_raster=miss))
    rows.sort(key=lambda r: (-r["pop_sum"], r["tile"]))
    best = rows[0]; d = load(best["key"], best["tile"]); pop = px.exposure(d["E"], d["N"])["pop"]
    ce, cn = wmedian(d["E"], pop), wmedian(d["N"], pop)
    bbox, (e0, n0) = tile_box_ll(best["tile"])
    out = dict(rule="max sum of GHSL 2020 POP at matched points over O' U N' (ties: lexicographic tile id)",
               tile=best["tile"], set=best["set"], n_points=best["n"], pop_sum=best["pop_sum"],
               window_centre_E=ce, window_centre_N=cn, window_side_m=WIN, tile_origin_E=e0, tile_origin_N=n0,
               osm_bbox_wgs84=bbox, n_tiles=len(rows), tiles_with_missing_raster=sum(r["n_no_raster"] > 0 for r in rows),
               top5=[{k: r[k] for k in ("tile", "set", "pop_sum", "n_no_raster")} for r in rows[:5]])
    os.makedirs(RES, exist_ok=True)
    json.dump(out, open(f"{RES}/selection.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


def osm_layers(tile):
    import pyproj
    from shapely.geometry import LineString
    js = json.load(open(f"{OSM}/{tile}_overpass.json", encoding="utf-8"))
    tr = pyproj.Transformer.from_crs(4326, 3035, always_xy=True)
    rail, mway, cities = [], [], []
    for el in js["elements"]:
        tg = el.get("tags", {})
        if el["type"] == "way" and "geometry" in el and len(el["geometry"]) >= 2:
            x, y = tr.transform([p["lon"] for p in el["geometry"]], [p["lat"] for p in el["geometry"]])
            ln = LineString(list(zip(x, y)))
            if tg.get("railway") == "rail" and tg.get("usage") == "main": rail.append(ln)
            elif tg.get("highway") == "motorway": mway.append(ln)
        elif el["type"] == "node" and tg.get("place") == "city":
            x, y = tr.transform(el["lon"], el["lat"])
            cities.append(dict(name=tg.get("name:en") or tg.get("name"), name_local=tg.get("name"), E=x, N=y))
    return rail, mway, cities


def near_lines(E, N, lines):
    from shapely import STRtree, points
    if not lines: return np.zeros(len(E), bool)
    tree = STRtree(lines); pts = points(np.c_[E, N])
    idx = tree.query(pts, predicate="dwithin", distance=BUF)
    m = np.zeros(len(E), bool); m[np.unique(idx[0])] = True
    return m


def states(a, b, T):
    A = np.abs(a) >= T; B = np.abs(b) >= T
    return np.select([A & B, A & ~B, ~A & B], [0, 1, 2], 3)   # 0 both, 1 A only (disappeared), 2 B only (new), 3 neither


def counts(d, m, T):
    dg = (d["gB"] - d["gA"])[m]
    s0 = states(d["vA"][m], d["vB"][m], T); s1 = states(d["vA"][m], d["vB"][m] - dg, T)
    f0 = np.isin(s0, (1, 2)); f1 = np.isin(s1, (1, 2))
    r = dict(n=int(m.sum()), dg_median=float(np.median(dg)) if m.any() else np.nan,
             dg_min=float(dg.min()) if m.any() else np.nan, dg_max=float(dg.max()) if m.any() else np.nan)
    for tag, s, f in (("C0", s0, f0), ("C1", s1, f1)):
        for k, lab in enumerate(("both", "A_only", "B_only", "neither")): r[f"{tag}_{lab}"] = int((s == k).sum())
        r[f"{tag}_flips"] = int(f.sum()); anyal = int((s != 3).sum())
        r[f"{tag}_alarm_any"] = anyal; r[f"{tag}_flip_rate"] = f.sum() / anyal if anyal else np.nan
    r["R"] = r["C1_flips"] / r["C0_flips"] if r["C0_flips"] else np.nan
    r["explained_by_dg"] = int((f0 & ~f1).sum()); r["new_after_dg"] = int((~f0 & f1).sum())
    return r, s0, s1


def stage_analyze():
    import pandas as pd
    sel = json.load(open(f"{RES}/selection.json")); t = sel["tile"]; key = [k for k, v in SETS.items() if v == sel["set"]][0]
    d = load(key, t); rail, mway, cities = osm_layers(t)
    ce, cn, h = sel["window_centre_E"], sel["window_centre_N"], WIN / 2
    areas = {"city_window": (np.abs(d["E"] - ce) <= h) & (np.abs(d["N"] - cn) <= h),
             "rail_main": near_lines(d["E"], d["N"], rail), "motorway": near_lines(d["E"], d["N"], mway),
             "tile": np.ones(len(d["E"]), bool)}
    city = min(cities, key=lambda c: (c["E"] - ce) ** 2 + (c["N"] - cn) ** 2) if cities else None
    rows, st = [], {}
    for T in TS:
        for a, m in areas.items():
            r, s0, s1 = counts(d, m, T); rows.append(dict(area=a, T=T, **r)); st[(a, T)] = (m, s0, s1)
    df = pd.DataFrame(rows); df.to_csv(f"{RES}/example_counts.csv", index=False)
    # representativeness (contract 4): percentile of the tile's T5 C0 pooled loss among O' U N' tiles (Confirm03 A_tile)
    A = pd.read_csv(f"{C03}/A_tile.csv")
    A = A[A.set.isin(["O'", "N'"]) & (A["var"] == "pub") & (A["mode"] == "C0") & (A["T"] == 5) & (A["dir"] == "loss") & (A.n_mov > 0)]
    A = A.assign(rate=A.n_lost / A.n_mov); me = A[A.tile == t]
    pct = float((A.rate < me.rate.iloc[0]).mean() * 100 + (A.rate == me.rate.iloc[0]).mean() * 50) if len(me) else np.nan
    whole = df[(df.area == "tile") & (df["T"] == 5)].iloc[0]
    summ = dict(tile=t, set=sel["set"], city=city, n_rail_ways=len(rail), n_motorway_ways=len(mway),
                percentile_T5_C0_loss=pct, n_tiles_percentile=int(len(A)),
                check_tile_loss_matches_A_tile=bool(len(me) and int(me.n_lost.iloc[0]) == int(whole["C0_A_only"])
                                                    and int(me.n_mov.iloc[0]) == int(whole["C0_both"] + whole["C0_A_only"])),
                counts=rows)
    json.dump(summ, open(f"{RES}/summary.json", "w"), indent=1, default=float)
    for T in TS: figure(d, st, T, t, city, rail, mway, (ce, cn), sel)
    print(df.to_string())
    print(json.dumps({k: v for k, v in summ.items() if k != "counts"}, indent=1, default=float))


def figure(d, st, T, t, city, rail, mway, ctr, sel):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    col = {0: "#2a78d6", 1: "#eb6834", 2: "#1baf7a", 3: "#c9c9c9"}
    lab = {0: "alarm in both", 1: "alarm only in 2020–2024", 2: "alarm only in 2020–2025", 3: "no alarm"}
    e0, n0 = sel["tile_origin_E"], sel["tile_origin_N"]; h = WIN / 2
    cname = city["name"] if city else "city window"
    titles = {"city_window": f"{cname}: 20 km window", "rail_main": "Main railways (100 m buffer)", "motorway": "Motorways (100 m buffer)"}
    fig, axs = plt.subplots(2, 3, figsize=(13, 8.6), constrained_layout=True)
    for j, a in enumerate(("city_window", "rail_main", "motorway")):
        m, s0, s1 = st[(a, T)]; E, N = d["E"][m] / 1000, d["N"][m] / 1000
        for i, (s, tag) in enumerate(((s0, "C0: as published"), (s1, "C1: 2020–2025 minus Δg"))):
            ax = axs[i, j]
            for ln in (rail if a != "motorway" else []) + (mway if a != "rail_main" else []):
                x, y = ln.xy; ax.plot(np.asarray(x) / 1000, np.asarray(y) / 1000, color="#555555", lw=0.6, zorder=1)
            for k, sz in ((3, 1.0), (0, 3.0), (1, 6.0), (2, 6.0)):
                q = s == k; ax.scatter(E[q], N[q], s=sz, c=col[k], lw=0, zorder=2 + (k != 3))
            if a == "city_window": ax.set_xlim((ctr[0] - h) / 1000, (ctr[0] + h) / 1000); ax.set_ylim((ctr[1] - h) / 1000, (ctr[1] + h) / 1000)
            else: ax.set_xlim(e0 / 1000, e0 / 1000 + 100); ax.set_ylim(n0 / 1000, n0 / 1000 + 100)
            ax.set_aspect("equal"); ax.tick_params(labelsize=7)
            nf = int(np.isin(s, (1, 2)).sum())
            ax.text(0.02, 0.98, f"flips: {nf:,}", transform=ax.transAxes, va="top", fontsize=8,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5))
            ax.set_title(f"{titles[a]}\n{tag}", fontsize=9)
            if j == 0: ax.set_ylabel("Northing (km, EPSG:3035)", fontsize=8)
            if i == 1: ax.set_xlabel("Easting (km, EPSG:3035)", fontsize=8)
    hd = [Line2D([], [], marker="o", ls="", color=col[k], label=lab[k], markersize=5) for k in (0, 1, 2, 3)]
    hd.append(Line2D([], [], color="#555555", lw=1, label="OSM lines"))
    fig.legend(handles=hd, loc="lower center", ncol=5, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.04))
    fig.suptitle(f"Tile {t} ({sel['set']}), vertical, |v| ≥ {T:g} mm/yr. Map data © OpenStreetMap contributors", fontsize=10)
    fig.savefig(f"{RES}/Fig_example_T{T:g}.png", dpi=200, bbox_inches="tight"); plt.close(fig)


def stage_smoke():
    for s, t in tiles()[:1] + [x for x in tiles() if x[0] == "Nprime"][:1]:
        d = load(s, t); print(s, t, {k: v.shape for k, v in d.items()})
    import p2lib_exposure as px; print("exposure fn:", px.exposure.__name__)
    for f in glob.glob(f"{OSM}/*_overpass.json"): print("osm file present:", os.path.basename(f))


if __name__ == "__main__":
    {"select": stage_select, "analyze": stage_analyze, "smoke": stage_smoke}[sys.argv[1]]()
