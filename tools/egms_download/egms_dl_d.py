"""EGMS downloader for Confirm03_Proximity part D (Windows). Reuses egms_dl.py (Confirm01) for the key, token, search
and download functions unchanged; adds two modes (contract section 5):
  ew      : L3 Ortho East-West files, releases 2020-2024 and 2020-2025, for the 208 tiles of O' and N'.
            Sub-release: the same sub_release as the U file used for that tile/release (manifests); else the highest.
  sparse  : set S'': continue the frozen Nt draw order of Confirm01 (N_draw_order.csv) after the last Nt rank used by N';
            skip tiles already in O or N' and tiles missing from the archive; download both U releases; compute
            n = min(nA, nB) epochs inside the common calendar window from the CSV headers; stop at 15 tiles with n <= 130
            or 40 tiles taken.
  probe_ew: find the API productType of the East-West product (tries known names; logged to probe_ew_log.txt).
The key is read from the environment variable as before (EGMS_SERVICE_KEY_JSON); its content, tokens and query ids are never
printed or written. Resumable: files already present are reused.
Usage: python egms_dl_d.py [probe_ew|ew|sparse|all]   (v6: Range resume of partial files; v5: parallel EW downloads; all = sparse first, then ew for N-prime then O-prime; writes DONE_D.txt when complete)"""
import csv, os, re, sys, zipfile, io
import egms_dl as base

HERE = base.HERE; PAPER = base.PAPER; OUTDIR = base.OUTDIR
DATA1 = os.path.join(PAPER, "contracts", "Confirm01_JAG", "data")
MAN = {"O'": os.path.join(DATA1, "manifest_Oprime.csv"), "N'": os.path.join(DATA1, "manifest_Nprime.csv")}
MAN_O = os.path.join(DATA1, "manifest_O.csv"); ORDER = os.path.join(DATA1, "N_draw_order.csv")
LOG = os.path.join(HERE, "download_log_D.csv"); PLOG = os.path.join(HERE, "probe_ew_log.txt")
RA, RB = "2020-2024", "2020-2025"; SPARSE_N, SPARSE_TARGET, SPARSE_MAX = 130, 15, 40
EW_CANDIDATES = ["ORTHO-EW", "ORTHO-EAST", "ORTHO-E", "ORTHO-EAST-WEST", "ORTHO-HORIZONTAL"]

class RobustSession(base.Session):
    """Token refresh with retries (the 12:37 run stopped on a read timeout at the token endpoint)."""
    def token(self):
        import time
        for attempt in range(6):
            try:
                with TOKLOCK: return super().token()
            except SystemExit: raise
            except Exception as e:
                w = 30 * (attempt + 1); print(f"  换取令牌时网络错误（{type(e).__name__}），{w} 秒后重试", flush=True); self.tok = None; time.sleep(w)
        sys.exit("多次无法换取令牌，已停止。重跑会续传。")

def search(sess, tile, rel, ptype, tolerant=False):
    import time
    for attempt in range(3):
        try:
            js = sess.call("search", {"id": None, "tileId": tile, "releases": [rel], "levels": ["L3"], "productType": ptype}); break
        except SystemExit:
            if tolerant: return None, []
            if attempt == 2: raise
            print(f"  搜索 {tile} {rel} 失败，60 秒后重试", flush=True); sess.tok = None; time.sleep(60)
    if not isinstance(js, dict) or not js.get("status"):
        return None, []
    return js.get("id"), [h for h in js.get("hits", []) if h.get("productType", ptype) == ptype]

def download(sess, tile, rel, fn, ptype, attempts=6):
    """v6: resumes partial files with an HTTP Range request (large O' files were failing mid-transfer and restarting
    from zero). If the server answers 206 the bytes are appended; if it answers 200 the file restarts from the start.
    The .part file is kept between attempts and between passes. Bearer token and fresh query id as in v2-v5."""
    import time, urllib.request, urllib.parse, urllib.error
    dst = os.path.join(OUTDIR, fn); part = dst + ".part"
    for attempt in range(attempts):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have and base.zip_ok(part): os.replace(part, dst); return "ok", os.path.getsize(dst), base.sha256(dst)
        qid, hits = search(sess, tile, rel, ptype)
        if not any(h["filename"] == fn for h in hits): return "vanished", 0, ""
        u = f"{base.API}/download/{urllib.parse.quote(fn)}?id={urllib.parse.quote(str(qid), safe='')}"
        hdr = {**base.UA, "Authorization": f"Bearer {sess.token()}"}
        if have: hdr["Range"] = f"bytes={have}-"
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers=hdr), timeout=180) as r:
                mode = "ab" if (have and r.status == 206) else "wb"
                if have: print(f"    {fn}: 从 {have/1e6:.0f} MB 处" + ("续传" if mode == "ab" else "重新开始（服务器不支持续传）"), flush=True)
                with open(part, mode) as f:
                    while True:
                        b = r.read(1 << 20)
                        if not b: break
                        f.write(b)
            if base.zip_ok(part): os.replace(part, dst); return "ok", os.path.getsize(dst), base.sha256(dst)
            print(f"    {fn}: 文件尚不完整，继续", flush=True)
        except urllib.error.HTTPError as e:
            import time as _t
            if e.code == 416 and base.zip_ok(part): os.replace(part, dst); return "ok", os.path.getsize(dst), base.sha256(dst)
            if e.code == 416: os.remove(part)
            if e.code == 401 and _t.time() - sess.t0 > 30 * 60: sess.tok = None
            w = 60 * (attempt + 1); print(f"    {fn}: HTTP {e.code}，{w} 秒后重试（第 {attempt + 1}/{attempts} 次）", flush=True); time.sleep(w)
        except Exception as e:
            got = os.path.getsize(part) if os.path.exists(part) else 0
            w = 30 if got > have else 60 * (attempt + 1)
            print(f"    {fn}: 网络中断（{type(e).__name__}），已有 {got/1e6:.0f} MB，{w} 秒后续传", flush=True); time.sleep(w)
    return "failed", 0, ""

def probe_ew(sess):
    lines = [f"{base.now()} 东西向产品类型探测（瓦片 {base.PROBE_TILE}，{RA}）"]; found = None
    for p in EW_CANDIDATES:
        _, hits = search(sess, base.PROBE_TILE, RA, p, tolerant=True)
        names = [h["filename"] for h in hits]
        lines.append(f"  {p}: " + (", ".join(names) if names else "无"))
        if found is None and any("_E_" in n for n in names): found = p
    lines.append(f"选用：{found}")
    open(PLOG, "w", encoding="utf-8").write("\n".join(lines) + "\n"); print("\n".join(lines), flush=True)
    if not found: sys.exit("没有找到东西向产品类型。请把 probe_ew_log.txt 发给 Claude。")
    return found

def sub_of(path):
    m = re.search(r"_(\d+)\.zip$", os.path.basename(path or "")); return int(m.group(1)) if m else None

def pick(hits, want_sub):
    if not hits: return None
    same = [h for h in hits if want_sub is not None and int(h.get("sub_release") or 0) == want_sub]
    return same[0] if same else base.best(hits)

class Log:
    def __init__(self):
        new = not os.path.exists(LOG); self.f = open(LOG, "a", newline="", encoding="utf-8"); self.w = csv.writer(self.f)
        if new: self.w.writerow(["utc", "part", "set_or_rank", "tile", "release", "file", "action", "bytes", "sha256", "path", "note"])
    def row(self, *a):
        with LOCK: self.w.writerow([base.now(), *a]); self.f.flush()

import threading
LOCK = threading.Lock(); TOKLOCK = threading.Lock()
WORKERS = int(os.environ.get("EGMS_WORKERS", "3"))
FAILS = [0]
def get(sess, idx, log, part, tag, t, rel, hit, ptype, note=""):
    try:
        p = _get(sess, idx, log, part, tag, t, rel, hit, ptype, note)
    except (SystemExit, Exception) as e:
        log.row(part, tag, t, rel, hit.get("filename", ""), "error", "", "", "", f"{type(e).__name__}"); p = None
        print(f"    {hit.get('filename','')}: 出错（{type(e).__name__}），跳过，重跑会续传", flush=True)
    with LOCK: FAILS[0] = 0 if p else FAILS[0] + 1
    if FAILS[0] >= 12: sys.exit("连续 12 个文件未成功，已停止。请把窗口最后几十行发给 Claude。")
    return p

def _get(sess, idx, log, part, tag, t, rel, hit, ptype, note=""):
    fn = hit["filename"]
    if fn in idx and base.zip_ok(idx[fn]):
        log.row(part, tag, t, rel, fn, "use_local", os.path.getsize(idx[fn]), "", idx[fn], note); return idx[fn]
    res, nb, h = download(sess, t, rel, fn, ptype)
    log.row(part, tag, t, rel, fn, "downloaded" if res == "ok" else "download_" + res, nb, h, os.path.join(OUTDIR, fn) if res == "ok" else "", note)
    if res == "ok": idx[fn] = os.path.join(OUTDIR, fn); return idx[fn]
    print(f"    {fn}: 未完成（{res}），重跑会续传", flush=True); return None

def run_ew(sess, ptype, sets=("N'", "O'")):
    """v5: N' first, then O'; files already on disk are skipped without an API call; WORKERS parallel downloads
    (default 3, env EGMS_WORKERS); each job = fresh search + download of one file."""
    from concurrent.futures import ThreadPoolExecutor
    os.makedirs(OUTDIR, exist_ok=True); idx = base.local_index(); log = Log(); jobs = []
    for name in sets:
        for r in csv.DictReader(open(MAN[name], encoding="utf-8")):
            t = r["tile"]
            for rel, upath in ((RA, r["pathA"]), (RB, r["pathB"])):
                want = sub_of(upath); guess = f"EGMS_L3_{t}_100km_E_{rel.replace('-', '_')}_{want}.zip"
                if guess in idx and base.zip_ok(idx[guess]): continue
                jobs.append((name, t, rel, want))
    print(f"东西向待下载 {len(jobs)} 个文件，并行 {WORKERS} 个。", flush=True)
    def job(a):
        name, t, rel, want = a
        try:
            _, hits = search(sess, t, rel, ptype)
        except SystemExit:
            log.row("D1", name, t, rel, "", "error", "", "", "", "search"); return False
        h = pick(hits, want)
        if h is None:
            log.row("D1", name, t, rel, "", "missing_in_archive", "", "", "", ""); print(f"[EW {name}] {t} {rel}: 存档中没有", flush=True); return True
        note = "same_sub_as_U" if want is not None and int(h.get("sub_release") or 0) == want else f"U_sub={want}"
        print(f"[EW {name}] {t} {rel}: {h['filename']}", flush=True)
        return get(sess, idx, log, "D1", name, t, rel, h, ptype, note) is not None
    with ThreadPoolExecutor(WORKERS) as ex: res = list(ex.map(job, jobs))
    missing = sum(1 for x in res if not x)
    print(f"东西向本轮结束：仍缺 {missing} 个文件。", flush=True); return missing

def header_dates(path):
    with zipfile.ZipFile(path) as z:
        n = [x for x in z.namelist() if x.endswith(".csv")][0]
        with z.open(n) as f: head = io.TextIOWrapper(f, encoding="utf-8").readline()
    return sorted(c for c in head.strip().split(",") if re.fullmatch(r"\d{8}", c))

def common_n(pa, pb):
    a, b = header_dates(pa), header_dates(pb)
    if not a or not b: return None
    lo, hi = max(a[0], b[0]), min(a[-1], b[-1])
    return min(sum(lo <= d <= hi for d in a), sum(lo <= d <= hi for d in b))

def run_sparse(sess):
    os.makedirs(OUTDIR, exist_ok=True); idx = base.local_index(); log = Log()
    order = sorted([r for r in csv.DictReader(open(ORDER, encoding="utf-8")) if r["stratum"] == "Nt"], key=lambda r: int(r["rank"]))
    nman = list(csv.DictReader(open(MAN["N'"], encoding="utf-8"))); used = {r["tile"] for r in nman}
    used_nt = {r["tile"] for r in nman if r.get("stratum") == "Nt"}
    last = max(int(r["rank"]) for r in order if r["tile"] in used_nt)
    skip = used | {r["tile"] for r in csv.DictReader(open(MAN_O, encoding="utf-8"))}
    taken = sparse = 0
    for r in order:
        if int(r["rank"]) <= last: continue
        if sparse >= SPARSE_TARGET or taken >= SPARSE_MAX: break
        t = r["tile"]
        if t in skip: log.row("D2", r["rank"], t, "", "", "skip_already_used", "", "", "", ""); continue
        local = {}
        for rel in (RA, RB):
            pre = f"EGMS_L3_{t}_100km_U_{rel.replace('-', '_')}_"
            c = sorted(f for f in idx if f.startswith(pre) and base.zip_ok(idx[f]))
            if c: local[rel] = idx[c[-1]]
        if len(local) == 2:
            paths = local
        else:
            plan = {}
            try:
                for rel in (RA, RB): _, hits = search(sess, t, rel, "ORTHO-UP"); plan[rel] = base.best(hits)
            except SystemExit:
                log.row("D2", r["rank"], t, "", "", "error", "", "", "", "search"); print(f"[S'' #{r['rank']}] {t}: 搜索失败，下一轮再试", flush=True); continue
            miss = [rel for rel in (RA, RB) if plan[rel] is None]
            if miss:
                log.row("D2", r["rank"], t, "/".join(miss), "", "skip_not_in_archive", "", "", "", ""); print(f"[S'' #{r['rank']}] {t}: 存档中缺 {'/'.join(miss)}，跳过", flush=True); continue
            paths = {rel: get(sess, idx, log, "D2", r["rank"], t, rel, plan[rel], "ORTHO-UP") for rel in (RA, RB)}
        if not all(paths.values()): print(f"[S'' #{r['rank']}] {t}: 下载未完成，重跑会续传；本瓦片暂不计数", flush=True); continue
        n = common_n(paths[RA], paths[RB]); taken += 1; sparse += int(n is not None and n <= SPARSE_N)
        log.row("D2", r["rank"], t, "both", "", "taken", "", "", "", f"n_header={n};taken={taken};sparse={sparse}")
        print(f"[S'' #{r['rank']}] {t}: n={n}；已取 {taken}，其中稀疏 {sparse}/{SPARSE_TARGET}", flush=True)
    log.row("D2", "", "", "", "", "stop", "", "", "", f"taken={taken};sparse={sparse}")
    print(f"补充稀疏瓦片本轮结束：取 {taken} 个，其中稀疏 {sparse} 个。", flush=True)
    return sparse >= SPARSE_TARGET or taken >= SPARSE_MAX

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    var, key = base.load_key(); sess = RobustSession(key)
    print(f"{base.now()} 使用环境变量：{var}（内容不记录）", flush=True)
    done = os.path.join(HERE, "DONE_D.txt")
    ok_sparse, miss = True, 0
    if mode in ("sparse", "all"): ok_sparse = run_sparse(sess)
    if mode in ("probe_ew", "ew", "all"):
        ptype = "ORTHO-EAST" if os.path.exists(PLOG) and "选用：ORTHO-EAST" in open(PLOG, encoding="utf-8").read() else probe_ew(sess)
    if mode in ("ew", "all"): miss = run_ew(sess, ptype)
    if mode == "all" and ok_sparse and miss == 0:
        open(done, "w", encoding="utf-8").write(base.now() + " all part-D files present\n"); print("D 部分全部文件已齐。", flush=True)
