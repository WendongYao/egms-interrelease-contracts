"""Confirm02_Mechanism s93 (VM): M0 metadata audit (CONTRACT section 4, M0). For every tile of the O, O' and N' manifests of
Confirm01 and both releases: GNSS model version (<gnss><version> of the tile xml), production date, file-name suffix (_1/_2),
whether the CSV carries a gnss_velocity_u column. Reads only the xml and the CSV header; no velocities, no differences.
Usage: python s93_meta.py MNT_ROOT CONFIRM01_DATA_DIR OUT_CSV"""
import sys, re, zipfile, csv, io
ROOT, D1, OUT = sys.argv[1:4]
rows = []
for name, man in (("O", "manifest_O.csv"), ("O'", "manifest_Oprime.csv"), ("N'", "manifest_Nprime.csv")):
    for r in csv.DictReader(open(f"{D1}/{man}")):
        for side in ("A", "B"):
            p = r[f"path{side}"]; z = zipfile.ZipFile(f"{ROOT}/{p}")
            x = z.read([n for n in z.namelist() if n.endswith(".xml")][0]).decode()
            ver = re.search(r"<gnss>\s*<version>(.*?)</version>", x, re.S); dat = re.search(r"<production_date>(.*?)</production_date>", x)
            cn = [n for n in z.namelist() if n.endswith(".csv")][0]
            head = io.TextIOWrapper(z.open(cn)).readline().strip().split(",")
            m = re.search(r"U_(\d{4})_(\d{4})_(\d+)\.zip$", p)
            rows.append(dict(set=name, tile=r["tile"], side=side, release=f"{m.group(1)}-{m.group(2)}", suffix=m.group(3),
                             gnss_version=ver.group(1).strip() if ver else "", production_date=dat.group(1) if dat else "",
                             has_gnss_u=int("gnss_velocity_u" in head), file=p.split("/")[-1]))
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("rows", len(rows))
