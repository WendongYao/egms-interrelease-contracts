"""Release Stability Check (RSC) - command line (Confirm02_Mechanism, CONTRACT section 4 Part III).
Input: two EGMS Ortho (L3 U) tile zips or CSVs of the same area from two releases (A = older, B = newer); optional AOI box
(EPSG:3035 metres) and optional GHSL directory for exposure weights.
Steps: 1 metadata check (GNSS model version, production date, file suffix); 2 reference harmonization (C1 with the per-point
gnss_velocity_u columns when both files carry them, otherwise C1c with a user constant, default 0 = no harmonization and a
warning); 3 optional tile-plane harmonization (--plane; off by default); 4 stability report (A1-A3 and per-hotspot labels
both / only_A / only_B) as JSON + CSV. Published velocities (mean_velocity) only. The report never says which release is better.
Usage: python rsc.py A.zip B.zip OUT_PREFIX [--aoi E0 N0 E1 N1] [--ghsl DIR] [--const C] [--plane]"""
import sys, os, re, io, json, zipfile, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rsc_core as rc
def read(path):
    if path.endswith(".zip"):
        z = zipfile.ZipFile(path); cn = [n for n in z.namelist() if n.endswith(".csv")][0]; xn = [n for n in z.namelist() if n.endswith(".xml")]
        xml = z.read(xn[0]).decode() if xn else ""; head = pd.read_csv(z.open(cn), nrows=0).columns
        use = [c for c in ("easting", "northing", "mean_velocity", "gnss_velocity_u") if c in head]
        df = pd.read_csv(z.open(cn), usecols=use, dtype=np.float64)
    else:
        xml = ""; head = pd.read_csv(path, nrows=0).columns; use = [c for c in ("easting", "northing", "mean_velocity", "gnss_velocity_u") if c in head]
        df = pd.read_csv(path, usecols=use, dtype=np.float64)
    ver = re.search(r"<gnss>\s*<version>(.*?)</version>", xml, re.S); dat = re.search(r"<production_date>(.*?)</production_date>", xml)
    suf = re.search(r"_(\d+)\.(zip|csv)$", path)
    meta = dict(file=os.path.basename(path), gnss_version=ver.group(1).strip() if ver else None, production_date=dat.group(1) if dat else None,
                suffix=suf.group(1) if suf else None, has_gnss_u="gnss_velocity_u" in df.columns)
    return df, meta
def run(pa, pb, out, aoi=None, ghsl=None, const=0.0, use_plane=False):
    A, ma = read(pa); B, mb = read(pb); notes = []
    if ma["gnss_version"] != mb["gnss_version"]:
        notes.append(f"GNSS model version differs ({ma['gnss_version']} -> {mb['gnss_version']}): the reference frame or model changed between releases "
                     "(EGMS GNSS Calibration Report v5.2, Table 3). Compare after harmonization.")
    ia, ib = rc.align(A.easting.values, A.northing.values, B.easting.values, B.northing.values)
    E = np.rint(A.easting.values[ia]).astype(np.int64); N = np.rint(A.northing.values[ia]).astype(np.int64)
    vA = A.mean_velocity.values[ia]; vB = B.mean_velocity.values[ib]
    if aoi is not None:
        m = (E >= aoi[0]) & (N >= aoi[1]) & (E < aoi[2]) & (N < aoi[3]); E, N, vA, vB = E[m], N[m], vA[m], vB[m]; ia, ib = ia[m], ib[m]
    if ma["has_gnss_u"] and mb["has_gnss_u"]:
        dg = B.gnss_velocity_u.values[ib] - A.gnss_velocity_u.values[ia]; mode = "C2" if use_plane else "C1"; vBh = rc.harmonize(E, N, vA, vB, dg, mode)
    else:
        dg = None; mode = "C1c"; vBh = rc.harmonize(E, N, vA, vB, mode="C1c", const=const)
        if use_plane: vBh = vBh - rc.plane(E, N, vBh - vA); mode = "C1c+plane"
        notes.append("gnss_velocity_u missing in at least one release: per-point harmonization impossible; constant used = %g mm/yr." % const)
    built = pop = None
    if ghsl:
        os.environ["GHSL_DIR"] = ghsl; import p2lib_exposure as px; w = px.exposure(E, N); built, pop = w["built"], w["pop"]
    res0, _ = rc.decision_metrics(E, N, vA, vB, built, pop, hotspot_table=True)
    res, tabs = rc.decision_metrics(E, N, vA, vBh, built, pop, hotspot_table=True)
    rep = dict(metadata=dict(A=ma, B=mb), n_points=int(len(E)), mode=mode, notes=notes,
               median_dg=float(np.median(dg)) if dg is not None else None, as_published=res0, harmonized=res)
    json.dump(rep, open(out + "_report.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    rows = [dict(threshold=th, **r) for th, rr in tabs.items() for r in rr]; pd.DataFrame(rows).to_csv(out + "_hotspots.csv", index=False)
    return rep, (E, N, vA, vB, dg)
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("A"); ap.add_argument("B"); ap.add_argument("out")
    ap.add_argument("--aoi", nargs=4, type=float); ap.add_argument("--ghsl"); ap.add_argument("--const", type=float, default=0.0); ap.add_argument("--plane", action="store_true")
    a = ap.parse_args(); rep, _ = run(a.A, a.B, a.out, a.aoi, a.ghsl, a.const, a.plane); print(json.dumps(dict(mode=rep["mode"], n=rep["n_points"], notes=rep["notes"]), indent=1))
