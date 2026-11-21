"""Copy completed EGMS east-west zips from the user's Downloads folder into raw_new (Confirm03 part D, browser route).
Only files listed in need_ew.json are copied (deleting from Downloads is not permitted here; originals stay there), after a full zip CRC test; others are left untouched. Prints what remains."""
import json, os, shutil, zipfile, hashlib, csv, sys
from datetime import datetime, timezone
HOME = os.path.expanduser("~"); DL = os.path.join(HOME, "mnt", "Downloads")
P = os.path.join(HOME, "mnt", "Paper_EGMSRelQuality_20261117"); RAW = os.path.join(P, "raw_new"); T = os.path.join(P, "tools", "egms_download")
need = json.load(open(os.path.join(T, "need_ew.json"))); want = {x["file"]: x for x in need}
log = os.path.join(T, "download_log_D.csv"); moved = []; MAXN = int(sys.argv[1]) if len(sys.argv) > 1 else 8
for f in sorted(os.listdir(DL)):
    if f not in want: continue
    src = os.path.join(DL, f); dst = os.path.join(RAW, f)
    if os.path.exists(dst): continue
    try:
        with zipfile.ZipFile(src) as z: bad = z.testzip()
    except Exception as e: print("not ready/invalid:", f, type(e).__name__); continue
    if bad: print("CRC error:", f, bad); continue
    if len(moved) >= MAXN: break
    n = os.path.getsize(src); tmp = dst + ".tmp"; shutil.copyfile(src, tmp)
    if os.path.getsize(tmp) != n: print("size mismatch", f); continue
    os.replace(tmp, dst); h = ""; moved.append(f)
    with open(log, "a", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerow([datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S"), "D1", want[f]["set"].replace("prime", "'"), want[f]["tile"], want[f]["rel"].replace("_", "-"), f, "downloaded_browser", n, h, dst, "Explorer"])
have = set(os.listdir(RAW)); rem = [x for x in need if x["file"] not in have]
print("moved", len(moved), moved); print("remaining", len(rem), "N'", sum(x["set"] == "Nprime" for x in rem))
