"""Confirm01_JAG s85: manifest of the O set (136 included RelCons01 tiles), file choice identical to RelCons01 s52
(sorted recursive glob of EGMS_L3_*_100km_U_20*_1.zip under raw_cache; first hit per tile and release).
Paths are written relative to the VM mount root (~/mnt), e.g. raw_cache/...
Usage: python s85_manifest_O.py RAW_CACHE BY_TILE_CSV OUT_CSV"""
import glob, re, sys, pandas as pd
root, bt, out = sys.argv[1:4]
by = {}
for f in sorted(glob.glob(f"{root}/**/EGMS_L3_*_100km_U_20*_1.zip", recursive=True)):
    m = re.search(r"EGMS_L3_(E\d+N\d+)_100km_U_(2019_2023|2020_2024)_1\.zip", f)
    if m: by.setdefault(m.group(1), {}).setdefault(m.group(2), f)
inc = pd.read_csv(bt); inc = inc[inc.included == True].tile.tolist(); assert len(inc) == 136
import os
rel = lambda p: "raw_cache/" + os.path.relpath(p, root).replace(os.sep, "/")
pd.DataFrame([dict(tile=t, set="O", pathA=rel(by[t]["2019_2023"]), pathB=rel(by[t]["2020_2024"])) for t in sorted(inc)]).to_csv(out, index=False)
print(len(inc))
