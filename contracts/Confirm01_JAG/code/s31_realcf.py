"""Probe_RealCF01 — truth change between EGMS versions under EGMSQA's deterministic rules (contract §3–5).
Rules ported verbatim from EGMSQA egms_autolabel_trend_qa_v5_8.py (year_fracs, fit_linear, fit_harmonic,
trend_dir_label, seasonality_label, detect_step)."""
import sys, json, math
from datetime import datetime
import numpy as np, pandas as pd

W = sys.argv[1]; R = f"{W}/results"
TOL_V, TOL_AMP, TOL_STEP = 0.2, 0.5, 0.5
DATE_TOL_DAYS = 7

# ---- ported EGMSQA rules ----
def year_fracs(date_cols):
    dts = [datetime.strptime(d, "%Y%m%d") for d in date_cols]; t0 = dts[0]
    return np.array([(dt - t0).days / 365.25 for dt in dts], dtype=np.float64)
def fit_linear(t, y):
    A = np.vstack([t, np.ones_like(t)]).T; coef, *_ = np.linalg.lstsq(A, y, rcond=None); return float(coef[0]), float(coef[1])
def fit_harmonic(t, y, period_years=1.0):
    w = 2.0 * math.pi / period_years
    A = np.vstack([np.sin(w * t), np.cos(w * t), np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None); a, b = float(coef[0]), float(coef[1]); return math.sqrt(a * a + b * b)
def trend_dir_label(s, tol=TOL_V):
    return "uplifting" if s >= tol else ("subsiding" if s <= -tol else "stable")
def seasonality_label(amp, tol=TOL_AMP):
    return "seasonal" if amp >= tol else "nonseasonal"
def detect_step(y, date_cols, rmse_mm, min_jump_mm=TOL_STEP):
    y = np.asarray(y, dtype=np.float64); finite = np.isfinite(y)
    if len(y) < 2 or finite.sum() < 2:
        return "no_step", "none", 0.0
    dy = np.diff(y); ok = finite[:-1] & finite[1:]
    if ok.sum() == 0:
        return "no_step", "none", 0.0
    thr = max(float(min_jump_mm), 3.0 * float(rmse_mm)); idxs = np.where(ok)[0]
    best = int(idxs[np.argmax(np.abs(dy[idxs]))]); mag = float(dy[best])
    if abs(mag) < thr or best + 1 >= len(date_cols):
        return "no_step", "none", 0.0
    return "step", str(date_cols[best + 1]), mag

def rmse_lin_harm(t, y):
    A = np.vstack([t, np.sin(2 * math.pi * t), np.cos(2 * math.pi * t), np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None); r = y - A @ coef; return float(np.sqrt(np.mean(r ** 2)))

def labels(y, dates, rmse_mm):
    y = np.asarray(y, dtype=np.float64); t = year_fracs(list(dates)); fin = np.isfinite(y)
    if fin.sum() < 3:
        return None
    v, _ = fit_linear(t[fin], y[fin]); amp = fit_harmonic(t[fin], y[fin])
    sb, sd, sm = detect_step(y, list(dates), rmse_mm)
    return dict(v=v, trend=trend_dir_label(v), amp=amp, seas=seasonality_label(amp), step=sb, step_date=sd, step_mag=sm, n=int(fin.sum()))

def days_between(d1, d2):
    return abs((datetime.strptime(d1, "%Y%m%d") - datetime.strptime(d2, "%Y%m%d")).days)

def compare(LA, LB):
    """Returns dict of 6 items: changed / substantive flags + magnitudes."""
    out = {}
    out["A_class"] = dict(changed=LA["trend"] != LB["trend"], subst=(LA["trend"] != LB["trend"]) and abs(LA["v"]) >= TOL_V and abs(LB["v"]) >= TOL_V,
                          sign_flip=(LA["trend"] != LB["trend"]) and LA["trend"] != "stable" and LB["trend"] != "stable", d=abs(LA["v"] - LB["v"]))
    dv = abs(LA["v"] - LB["v"]); out["A_num"] = dict(changed=dv > TOL_V, subst=dv > 2 * TOL_V, d=dv)
    da = abs(LA["amp"] - LB["amp"])
    out["B_class"] = dict(changed=LA["seas"] != LB["seas"], subst=(LA["seas"] != LB["seas"]) and da >= TOL_AMP, d=da)
    out["B_num"] = dict(changed=da > TOL_AMP, subst=da > 2 * TOL_AMP, d=da)
    sc = LA["step"] != LB["step"]
    appear_mag = abs(LB["step_mag"]) if LB["step"] == "step" else abs(LA["step_mag"])
    out["C_class"] = dict(changed=sc, subst=sc and appear_mag >= 1.0, d=appear_mag if sc else 0.0)
    if LA["step"] == "step" and LB["step"] == "step":
        dm = abs(LA["step_mag"] - LB["step_mag"]); dd = days_between(LA["step_date"], LB["step_date"])
        out["C_num"] = dict(changed=(dm > TOL_STEP) or (dd > DATE_TOL_DAYS), subst=(dm > 2 * TOL_STEP) or (dd > DATE_TOL_DAYS), d=dm, dd=dd)
    elif sc:
        out["C_num"] = dict(changed=True, subst=appear_mag >= 1.0, d=appear_mag, dd=np.nan)
    else:
        out["C_num"] = dict(changed=False, subst=False, d=0.0, dd=np.nan)
    return out

def main():
    d = np.load(f"{W}/data/e1_panel.npz", allow_pickle=True)
    tiles = sorted({k.split("__")[0] for k in d.files})
    rows = []; tile_info = []
    for tile in tiles:
        A, B, dA, dB = d[f"{tile}__A"], d[f"{tile}__B"], list(d[f"{tile}__dA"]), list(d[f"{tile}__dB"])
        lo, hi = max(dA[0], dB[0]), min(dA[-1], dB[-1])
        ia = [i for i, x in enumerate(dA) if lo <= x <= hi]; ib = [i for i, x in enumerate(dB) if lo <= x <= hi]
        tile_info.append(dict(tile=tile, n_points=A.shape[0], W_start=lo, W_end=hi, nA_epochs=len(ia), nB_epochs=len(ib)))
        if len(ia) < 3 or len(ib) < 3:
            continue
        dAw = [dA[i] for i in ia]; dBw = [dB[i] for i in ib]
        for p in range(A.shape[0]):
            ya, yb = A[p, ia].astype(float), B[p, ib].astype(float)
            fin = np.isfinite(ya)
            if fin.sum() < 3:
                continue
            rmse = rmse_lin_harm(year_fracs(dAw)[fin], ya[fin])  # version-A based, shared threshold
            for win in ["W", "W64"]:
                if win == "W64":
                    ya2, da2 = ya[-64:], dAw[-64:]; yb2, db2 = yb[-64:], dBw[-64:]
                else:
                    ya2, da2, yb2, db2 = ya, dAw, yb, dBw
                LA, LB = labels(ya2, da2, rmse), labels(yb2, db2, rmse)
                if LA is None or LB is None:
                    continue
                cmp_ = compare(LA, LB)
                for item, c in cmp_.items():
                    rows.append(dict(tile=tile, point=p, window=win, item=item, changed=bool(c["changed"]), subst=bool(c["subst"]),
                                     sign_flip=bool(c.get("sign_flip", False)), d=float(c["d"]), dd=c.get("dd", np.nan),
                                     vA=LA["v"], vB=LB["v"], ampA=LA["amp"], ampB=LB["amp"], stepA=LA["step"], stepB=LB["step"], rmseA=rmse))
    per = pd.DataFrame(rows); per.to_csv(f"{R}/per_item.csv", index=False)
    pd.DataFrame(tile_info).to_csv(f"{R}/tile_windows.csv", index=False)
    S = {}
    for win in ["W", "W64"]:
        sub = per[per.window == win]
        fam = sub.groupby("item").agg(n=("changed", "size"), changed=("changed", "sum"), subst=("subst", "sum"), sign_flip=("sign_flip", "sum"), med_d_changed=("d", lambda x: float(np.median(x[sub.loc[x.index, "changed"]])) if sub.loc[x.index, "changed"].any() else np.nan)).reset_index()
        fam["changed_rate"] = fam.changed / fam.n
        fam.to_csv(f"{R}/by_family_{win}.csv", index=False)
        n_items = len(sub); n_changed = int(sub.changed.sum()); n_subst = int(sub.subst.sum())
        sign_flips = int(sub[(sub.item == "A_class")].sign_flip.sum())
        step_app = int(sub[(sub.item == "C_class")].subst.sum())
        S[win] = dict(n_points=int(sub.groupby(["tile", "point"]).ngroups), n_items=n_items, n_changed=n_changed, changed_rate=n_changed / n_items if n_items else np.nan,
                      n_substantive=n_subst, sign_flips_subst=sign_flips, step_appear_disappear_subst=step_app,
                      by_family=fam.to_dict("records"))
    tb = per[per.window == "W"].groupby("tile").agg(n_items=("changed", "size"), changed_rate=("changed", "mean"), subst=("subst", "sum"), rmseA_med=("rmseA", "median")).reset_index()
    tb = tb.merge(pd.DataFrame(tile_info), on="tile"); tb.to_csv(f"{R}/by_tile.csv", index=False)
    Wd = S["W"]
    S["P1_rate_ge_5pct"] = bool(Wd["changed_rate"] >= 0.05); S["P2_substantive_ge_500"] = bool(Wd["n_substantive"] >= 500)
    S["P3_discrete_ge_100"] = bool(Wd["sign_flips_subst"] >= 100 or Wd["step_appear_disappear_subst"] >= 100)
    S["decision"] = "PASS -> LLM experiment design" if (S["P1_rate_ge_5pct"] and S["P2_substantive_ge_500"] and S["P3_discrete_ge_100"]) else "STOP"
    with open(f"{R}/summary.json", "w") as f:
        json.dump(S, f, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o), ensure_ascii=False)
    print(json.dumps({k: v for k, v in S.items() if k not in ("W", "W64")}, indent=2))
    for win in ["W", "W64"]:
        print(win, {k: v for k, v in S[win].items() if k != "by_family"})
        print(pd.DataFrame(S[win]["by_family"]).round(3).to_string())

if __name__ == "__main__":
    main()
