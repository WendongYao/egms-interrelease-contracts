"""Confirm05_Example: fetch OSM features for the selected tile once, from the user's VM (contract section 2).

Usage: python s100_osm_fetch.py <tile> <south> <west> <north> <east> [--dry-run]
Writes data/osm/<tile>_overpass.json, <tile>_query.overpassql and <tile>_meta.json (time, bytes, SHA-256).
Query is pinned to 2026-11-21T00:00:00Z so it can be repeated; only railway=rail ways, highway=motorway ways
and place=city nodes are requested.
"""
import sys, os, json, hashlib, time, urllib.request, urllib.parse, datetime

URL = "https://overpass-api.de/api/interpreter"
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(os.path.dirname(HERE), "data", "osm")


def query(s, w, n, e):
    b = f"{s:.6f},{w:.6f},{n:.6f},{e:.6f}"
    return ('[out:json][timeout:300][date:"2026-11-21T00:00:00Z"];\n'
            f'(\n  way["railway"="rail"]({b});\n  way["highway"="motorway"]({b});\n  node["place"="city"]({b});\n);\n'
            'out geom;\n')


def main():
    tile = sys.argv[1]; s, w, n, e = map(float, sys.argv[2:6]); q = query(s, w, n, e)
    os.makedirs(OUT, exist_ok=True)
    open(f"{OUT}/{tile}_query.overpassql", "w", encoding="utf-8").write(q)
    if "--dry-run" in sys.argv:
        print(q); return
    req = urllib.request.Request(URL, data=urllib.parse.urlencode({"data": q}).encode(),
                                 headers={"User-Agent": "EGMSRelQuality-Confirm05/1.0 (research; UCD)",
                                          "Accept": "application/json"})
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=330) as r:
        body = r.read(); status = r.status
    js = json.loads(body)
    open(f"{OUT}/{tile}_overpass.json", "wb").write(body)
    meta = dict(tile=tile, url=URL, status=status, bytes=len(body), sha256=hashlib.sha256(body).hexdigest(),
                fetched_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
                seconds=round(time.time() - t0, 1), n_elements=len(js.get("elements", [])),
                osm_base=js.get("osm3s", {}).get("timestamp_osm_base"), remark=js.get("remark"),
                query_sha256=hashlib.sha256(q.encode()).hexdigest())
    json.dump(meta, open(f"{OUT}/{tile}_meta.json", "w"), indent=1)
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
