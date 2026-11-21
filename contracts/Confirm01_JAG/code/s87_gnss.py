"""Confirm01_JAG s87 (cloud) = Probe_RelCons01 s53 generalised: NGL IGS20 tenv3 for stations with >= 3 near points
(read from a local cache; downloaded only if missing), cut to the tile's window W, fit intercept + trend + annual +
semiannual + offsets at steps.txt epochs (one 5xMAD outlier pass). Station rules identical to s53.
Usage: python s87_gnss.py SAMPLE_DIR TENV3_CACHE STEPS_TXT WORKDIR_OUT   Writes WORKDIR_OUT/gnss_fit.csv."""
import sys, os, glob, json, math, time, urllib.request, concurrent.futures as cf
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rc_common import window, d2mjd, MJD0
from datetime import datetime
SD, CACHE, STEPS, W = sys.argv[1:5]; G = f"{W}/gnss_w"; os.makedirs(G, exist_ok=True); os.makedirs(W, exist_ok=True)
URL = "https://geodesy.unr.edu/gps_timeseries/IGS20/tenv3/IGS20/{sta}.tenv3"
def ymd_to_mjd(s): return datetime.strptime(s, "%y%b%d").date().toordinal() - MJD0
steps = {}
with open(STEPS, errors="ignore") as f:
    for line in f:
        p = line.split()
        if len(p) < 3:
            continue
        try:
            steps.setdefault(p[0], set()).add(ymd_to_mjd(p[1]))
        except Exception:
            continue
jobs = []
for jf in sorted(glob.glob(f"{SD}/*.json")):
    m = json.load(open(jf))
    if "excluded" in m:
        continue
    z = np.load(jf.replace(".json", ".npz"), allow_pickle=True)
    lo, hi, ia, ib = window(z["dA"], z["dB"])
    st = pd.Series(z["sta"]).value_counts()
    for s, n in st.items():
        if n >= 3:
            jobs.append((s, m["tile"], d2mjd(lo), d2mjd(hi), int(n)))
print("stations with >=3 near points:", len(jobs), flush=True)
def fetch_cut(job):
    sta, tile, m0, m1, n = job; out = f"{G}/{sta}.npz"
    if os.path.exists(out):
        return sta, "cached"
    cp = f"{CACHE}/{sta}.tenv3"
    for attempt in range(3):
        try:
            if os.path.exists(cp) and os.path.getsize(cp) > 0:
                txt = open(cp, encoding="ascii", errors="ignore").read(); break
            req = urllib.request.Request(URL.format(sta=sta), headers={"User-Agent": "confirm01-jag/1.0"})
            txt = urllib.request.urlopen(req, timeout=120).read().decode("ascii", "ignore"); break
        except Exception as e:
            err = type(e).__name__; time.sleep(3 * (attempt + 1))
    else:
        np.savez(out, mjd=np.zeros(0, int), up=np.zeros(0), status=f"error:{err}"); return sta, f"error:{err}"
    mj, up, seen = [], [], set()
    for line in txt.splitlines()[1:]:
        p = line.split()
        if len(p) < 17:
            continue
        try:
            mm = int(p[3]); u = (float(p[11]) + float(p[12])) * 1000.0
        except ValueError:
            continue
        if mm in seen or mm < m0 or mm > m1:
            continue
        seen.add(mm); mj.append(mm); up.append(u)
    o = np.argsort(mj); np.savez(out, mjd=np.array(mj, int)[o], up=np.array(up)[o], status="ok"); return sta, "ok"
with cf.ThreadPoolExecutor(6) as ex:
    st = list(ex.map(fetch_cut, jobs))
print("download status:", pd.Series([s for _, s in st]).value_counts().to_dict(), flush=True)
rows = []
for sta, tile, m0, m1, n in jobs:
    z = np.load(f"{G}/{sta}.npz", allow_pickle=True); mjd, up = z["mjd"], z["up"]
    span = m1 - m0 + 1; compl = len(mjd) / span
    sk = sorted(k for k in steps.get(sta, ()) if m0 < k <= m1)
    edge = any((k - m0) < 90 or (m1 - k) < 90 for k in sk)
    r = dict(sta=sta, tile=tile, n_near=n, n_days=int(len(mjd)), completeness=compl, n_steps=len(sk), step_near_edge=edge, status=str(z["status"]))
    ok = (compl >= 0.70) and (len(sk) <= 2) and (not edge) and len(mjd) > 30
    r["reason"] = "" if ok else ("download" if len(mjd) == 0 else ("completeness" if compl < 0.70 else ("steps>2" if len(sk) > 2 else "step_near_edge")))
    if ok:
        t = (mjd - d2mjd("20200101")) / 365.25
        def design(tt, mm):
            cols = [np.ones_like(tt), tt, np.sin(2 * math.pi * tt), np.cos(2 * math.pi * tt), np.sin(4 * math.pi * tt), np.cos(4 * math.pi * tt)]
            cols += [(mm >= k).astype(float) for k in sk]
            return np.vstack(cols).T
        X = design(t, mjd); coef, *_ = np.linalg.lstsq(X, up, rcond=None); res = up - X @ coef
        med = np.median(res); mad = np.median(np.abs(res - med)); keep = np.abs(res - med) <= 5 * mad if mad > 0 else np.ones(len(res), bool)
        X2 = design(t[keep], mjd[keep]); coef, *_ = np.linalg.lstsq(X2, up[keep], rcond=None); res2 = up[keep] - X2 @ coef
        s2 = float(res2 @ res2) / max(len(res2) - X2.shape[1], 1); cov = s2 * np.linalg.pinv(X2.T @ X2)
        r.update(v_gnss=float(coef[1]), v_gnss_se_white=float(math.sqrt(max(cov[1, 1], 0))), n_outliers=int((~keep).sum()), wrms=float(np.sqrt(np.mean(res2 ** 2))))
    r["included"] = bool(ok); rows.append(r)
g = pd.DataFrame(rows); g.to_csv(f"{W}/gnss_fit.csv", index=False)
print("gnss fits written:", len(g), "included:", int(g.included.sum()), "reasons:", g.reason.value_counts().to_dict())
