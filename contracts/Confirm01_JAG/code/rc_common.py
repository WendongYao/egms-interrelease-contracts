"""Probe_RelCons01 shared definitions (contract section 3)."""
import math
from datetime import date, datetime
import numpy as np
MJD0 = date(1858, 11, 17).toordinal()
T_REF = date(2020, 1, 1).toordinal()
SEED = 20261110
def d2ord(s): return datetime.strptime(str(s), "%Y%m%d").date().toordinal()
def d2mjd(s): return d2ord(s) - MJD0
def window(dA, dB):
    """Common calendar window W and the epoch indices of each release inside it."""
    dA = [str(x) for x in dA]; dB = [str(x) for x in dB]
    lo, hi = max(dA[0], dB[0]), min(dA[-1], dB[-1])
    ia = [i for i, x in enumerate(dA) if lo <= x <= hi]; ib = [i for i, x in enumerate(dB) if lo <= x <= hi]
    return lo, hi, ia, ib
def tyears(dates): return np.array([(d2ord(x) - T_REF) / 365.25 for x in dates], dtype=np.float64)
def fit_primary(t, y, min_n=60):
    """intercept + trend + annual sin/cos; returns (v, formal_se) or None."""
    fin = np.isfinite(y)
    if fin.sum() < min_n:
        return None
    t = t[fin]; y = y[fin].astype(np.float64)
    X = np.vstack([np.ones_like(t), t, np.sin(2 * math.pi * t), np.cos(2 * math.pi * t)]).T
    coef, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ coef
    s2 = float(r @ r) / max(len(y) - 4, 1); cov = s2 * np.linalg.pinv(X.T @ X)
    return float(coef[1]), float(math.sqrt(max(cov[1, 1], 0.0)))
def fit_ols(t, y, min_n=60):
    fin = np.isfinite(y)
    if fin.sum() < min_n:
        return None
    t = t[fin]; y = y[fin].astype(np.float64)
    X = np.vstack([np.ones_like(t), t]).T; coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return float(coef[1])
def density_class(n): return "dense" if n >= 200 else ("sparse" if n <= 130 else "mid")
