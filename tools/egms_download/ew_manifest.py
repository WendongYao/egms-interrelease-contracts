"""Confirm03 part D-1: locate every east-west (E) file needed for O' and N' and record its size and SHA-256.
Expected file per tile and release: same sub-release as the U file in the Confirm01 manifests (contract D-1).
Location: raw_new/ (API route) first, else ~/mnt/Downloads (browser route; files there are only read, never moved).
Outputs (resumable; rows already hashed are kept):
  tools/egms_download/ew_files_manifest.csv  set,tile,rel,file,location,bytes,sha256
  contracts/Confirm03_Proximity/data/manifest_E_{Oprime,Nprime}.csv  tile,set,pathA,pathB (relative to ~/mnt), as s82 expects
Usage: python3 ew_manifest.py [time_budget_s]
Prints counts only."""
import hashlib, os, re, sys, time
import pandas as pd

T0 = time.time(); BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 160
MNT = os.path.expanduser("~/mnt"); P = os.path.join(MNT, "Paper_EGMSRelQuality_20261117")
RAW = os.path.join(P, "raw_new"); DL = os.path.join(MNT, "Downloads")
OUTF = os.path.join(P, "tools", "egms_download", "ew_files_manifest.csv")
DD = os.path.join(P, "contracts", "Confirm03_Proximity", "data"); os.makedirs(DD, exist_ok=True)
sub = lambda p: re.search(r"_(\d+)\.zip$", p).group(1)

rows = []
for tag, setname in (("Oprime", "O'"), ("Nprime", "N'")):
    m = pd.read_csv(os.path.join(P, "contracts", "Confirm01_JAG", "data", f"manifest_{tag}.csv"))
    for _, r in m.iterrows():
        for rel, up in (("2020_2024", r.pathA), ("2020_2025", r.pathB)):
            rows.append(dict(set=setname, tag=tag, tile=r.tile, rel=rel, file=f"EGMS_L3_{r.tile}_100km_E_{rel}_{sub(up)}.zip"))
need = pd.DataFrame(rows)

done = pd.read_csv(OUTF) if os.path.exists(OUTF) else pd.DataFrame(columns=["set", "tile", "rel", "file", "location", "bytes", "sha256"])
have = set(done.file)
new = []; missing = 0; stopped = False
for _, r in need.iterrows():
    if r.file in have: continue
    if time.time() - T0 > BUDGET: stopped = True; break
    for loc, d in (("raw_new", RAW), ("Downloads", DL)):
        p = os.path.join(d, r.file)
        if os.path.exists(p): break
    else:
        missing += 1; continue
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    new.append(dict(set=r.set, tile=r.tile, rel=r.rel, file=r.file, location=loc, bytes=os.path.getsize(p), sha256=h.hexdigest()))
out = pd.concat([done, pd.DataFrame(new)], ignore_index=True) if new else done
out.to_csv(OUTF, index=False)
print(f"hashed now {len(new)}; total {len(out)}/{len(need)}; missing {missing}; stopped_on_budget {stopped}")

if len(out) == len(need):
    loc = dict(zip(out.file, out.location))
    rel_path = lambda f: (f"Paper_EGMSRelQuality_20261117/raw_new/{f}" if loc[f] == "raw_new" else f"Downloads/{f}")
    for tag in ("Oprime", "Nprime"):
        n = need[need.tag == tag]
        a = n[n.rel == "2020_2024"].set_index("tile"); b = n[n.rel == "2020_2025"].set_index("tile")
        man = pd.DataFrame(dict(tile=a.index, set=a.set.values, pathA=[rel_path(f) for f in a.file], pathB=[rel_path(b.loc[t, "file"]) for t in a.index]))
        man.to_csv(os.path.join(DD, f"manifest_E_{tag}.csv"), index=False)
        print(f"manifest_E_{tag}.csv: {len(man)} tiles")
