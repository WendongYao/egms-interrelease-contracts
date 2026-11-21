"""Confirm01_JAG s92 (cloud) = Probe_RelCons01 s51 generalised to a tile list: NGL stations (DataHoldings.txt,
region box 25<=lat<=75, -35<=lon<=50) -> EPSG:3035 -> stations inside the given tiles. Tile naming as in EGMS
file names (E<e>N<n>, no zero padding). Usage: python s92_stations.py DataHoldings.txt TILES.txt OUT.csv"""
import sys, numpy as np, pandas as pd, pyproj
DH, TL, OUT = sys.argv[1:4]
tiles = set(x.strip() for x in open(TL) if x.strip())
rows = []
with open(DH, errors="ignore") as f:
    f.readline()
    for line in f:
        p = line.split()
        if len(p) < 9: continue
        try: lat, lon, h = float(p[1]), float(p[2]), float(p[3])
        except ValueError: continue
        if lon > 180: lon -= 360
        if 25 <= lat <= 75 and -35 <= lon <= 50: rows.append((p[0], lat, lon, h, p[7], p[8]))
s = pd.DataFrame(rows, columns=["sta", "lat", "lon", "h", "dtbeg", "dtend"]).drop_duplicates("sta")
E, N = pyproj.Transformer.from_crs(4326, 3035, always_xy=True).transform(s.lon.values, s.lat.values); s["E"] = E; s["N"] = N
s["tile"] = ["E%dN%d" % (int(e // 1e5), int(n // 1e5)) for e, n in zip(E, N)]
s = s[s.tile.isin(tiles)].sort_values(["tile", "sta"]).reset_index(drop=True)
s.to_csv(OUT, index=False); print("stations:", len(s), "tiles with >=1 station:", s.tile.nunique(), "of", len(tiles))
