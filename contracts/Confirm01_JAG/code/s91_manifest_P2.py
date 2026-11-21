"""Confirm01_JAG s91: manifests for the P2 pair (Amendment 1): N' (new tiles) and O' (old tiles), from the
downloader's log (last state per stratum/tile/release). Paths relative to ~/mnt. Tiles with a missing or failed
file are listed in *_pending.csv and left out until the download is completed (the manifest is rebuilt then).
Usage: python s91_manifest_P2.py PAPER_DIR"""
import csv, os, sys, pandas as pd
P = sys.argv[1]; D = f"{P}/contracts/Confirm01_JAG/data"
rows = list(csv.DictReader(open(f"{P}/tools/egms_download/download_log.csv", encoding="utf-8")))
mo = pd.read_csv(f"{D}/manifest_O.csv").set_index("tile")
last = {}
for r in rows:
    if r["file"]: last[(r["stratum"], r["tile"], r["release"])] = r
def rel_path(r, tile):
    a = r["action"]
    if a == "downloaded": return "Paper_EGMSRelQuality_20261117/raw_new/" + r["file"]
    if a == "use_local":
        p = r["path"].replace("\\", "/"); k = p.find("raw_cache/")
        if k >= 0: return p[k:]
        k = p.find("Paper_EGMSRelQuality_20261117/raw_new/")
        return p[k:] if k >= 0 else None  # files fetched by an earlier run are logged as use_local on re-runs
    if a == "use_local_P1B": return mo.loc[tile, "pathB"]
    return None
out = {"N": [], "O": []}; pend = {"N": [], "O": []}
strata = {}
for (st, t, rel), r in last.items(): strata.setdefault((st, t), {})[rel] = r
for (st, t), d in sorted(strata.items()):
    g = "O" if st == "O'" else "N"
    pa = rel_path(d["2020-2024"], t) if "2020-2024" in d else None
    pb = rel_path(d["2020-2025"], t) if "2020-2025" in d else None
    if pa and pb: out[g].append(dict(tile=t, set=("O'" if g == "O" else "N'"), stratum=st, pathA=pa, pathB=pb))
    else: pend[g].append(dict(tile=t, stratum=st, has_A=bool(pa), has_B=bool(pb)))
for g, name in (("N", "Nprime"), ("O", "Oprime")):
    pd.DataFrame(out[g]).to_csv(f"{D}/manifest_{name}.csv", index=False)
    pd.DataFrame(pend[g]).to_csv(f"{D}/manifest_{name}_pending.csv", index=False)
    print(name, "complete", len(out[g]), "pending", len(pend[g]))
