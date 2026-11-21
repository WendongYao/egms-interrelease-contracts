"""EGMS Ortho-U downloader for Confirm01_JAG (Windows; needs either `pyjwt` or `cryptography` for RS256 signing).

Auth flow follows the official notebook github.com/copernicus-land/egms-api:
CLMS service key (JSON, read from an environment variable, never printed or written) -> signed JWT grant ->
short-lived access token at token_uri -> POST /insar-api/archive/search per tile -> download
/insar-api/archive/download/<filename>?id=<query id>. Query ids and tokens are never logged.

Mode `all` (default): probe first (both P2 releases must exist), then
  (1) N' set: N-set tiles in the frozen draw order with BOTH 2020-2024 and 2020-2025, until each quota is met;
  (2) O' set: for the 136 O tiles, 2020-2025 (and 2020-2024 only if a higher sub_release than the local file exists).
Resumable. (Amendment 1 to the Confirm01_JAG contract.)"""
import base64, csv, hashlib, json, os, re, sys, time, urllib.request, urllib.error, urllib.parse, zipfile
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))
ORDER = os.path.join(PAPER, "contracts", "Confirm01_JAG", "data", "N_draw_order.csv")
OUTDIR = os.path.join(PAPER, "raw_new")
RAW_CACHE = r"C:\Users\Jupyter\Desktop\DeformCLIP\data\raw_cache"
LOG = os.path.join(HERE, "download_log.csv")
PROBE_LOG = os.path.join(HERE, "probe_log.txt")
API = "https://egms.land.copernicus.eu/insar-api/archive"
QUOTA = {"Al": 12, "S": 20, "C": 20, "Nt": 20}
MAX_RANK = 45
PROBE_TILE = "E15N10"
RA, RB = "2020-2024", "2020-2025"  # Amendment 1: pair P2
MANIFEST_O = os.path.join(PAPER, "contracts", "Confirm01_JAG", "data", "manifest_O.csv")
UA = {"User-Agent": "confirm01-jag/2.0"}

def now(): return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
def b64u(b): return base64.urlsafe_b64encode(b).rstrip(b"=").decode()

def load_key():
    var = os.environ.get("EGMS_VAR") or "EGMS_SERVICE_KEY_JSON"
    raw = os.environ.get(var)
    if not raw: sys.exit(f"环境变量 {var} 不存在。")
    raw = raw.strip()
    if not raw.startswith("{") and os.path.exists(raw): raw = open(raw, encoding="utf-8").read()
    try: k = json.loads(raw)
    except Exception: sys.exit(f"环境变量 {var} 不是合法的 JSON（内容不显示）。")
    need = {"private_key", "client_id", "user_id", "token_uri"}
    if not need <= set(k): sys.exit(f"{var} 缺少字段：{sorted(need - set(k))}（只列字段名）。")
    return var, k

def sign_rs256(msg, pem):
    try:
        import jwt  # pyjwt
        return None  # handled in make_grant
    except ImportError:
        pass
    try:
        from cryptography.hazmat.primitives import hashes, serialization
        from cryptography.hazmat.primitives.asymmetric import padding
    except ImportError:
        sys.exit("当前 Python 环境既没有 pyjwt 也没有 cryptography，无法签名。请告诉 Claude，由你决定是否安装（pip install pyjwt cryptography）。")
    key = serialization.load_pem_private_key(pem.encode(), password=None)
    return key.sign(msg, padding.PKCS1v15(), hashes.SHA256())

def make_grant(k):
    t = int(time.time())
    claims = {"iss": k["client_id"], "sub": k["user_id"], "aud": k["token_uri"], "iat": t, "exp": t + 3600}
    try:
        import jwt
        g = jwt.encode(claims, k["private_key"].encode(), algorithm="RS256")
        return g if isinstance(g, str) else g.decode()
    except ImportError:
        head = b64u(json.dumps({"alg": "RS256", "typ": "JWT"}, separators=(",", ":")).encode())
        body = b64u(json.dumps(claims, separators=(",", ":")).encode())
        sig = sign_rs256(f"{head}.{body}".encode(), k["private_key"])
        return f"{head}.{body}.{b64u(sig)}"

class Session:
    def __init__(self, key): self.key = key; self.tok = None; self.t0 = 0
    def token(self):
        if self.tok and time.time() - self.t0 < 45 * 60: return self.tok
        data = urllib.parse.urlencode({"grant_type": "urn:ietf:params:oauth:grant-type:jwt-bearer", "assertion": make_grant(self.key)}).encode()
        req = urllib.request.Request(self.key["token_uri"], data=data, headers={**UA, "Accept": "application/json", "Content-Type": "application/x-www-form-urlencoded"})
        try:
            with urllib.request.urlopen(req, timeout=60) as r: js = json.loads(r.read())
        except urllib.error.HTTPError as e:
            sys.exit(f"换取访问令牌失败：HTTP {e.code}。密钥可能已失效，或需要在 CLMS 网站重新生成。")
        if not js.get("access_token"): sys.exit("换取访问令牌失败：返回中没有 access_token（内容不显示）。")
        self.tok = js["access_token"]; self.t0 = time.time(); return self.tok
    def call(self, path, body=None):
        for attempt in range(4):
            hdr = {**UA, "Authorization": f"Bearer {self.token()}", "Accept": "application/json"}
            if body is not None: hdr["Content-Type"] = "application/json"
            req = urllib.request.Request(f"{API}/{path}", data=None if body is None else json.dumps(body).encode(), headers=hdr)
            try:
                with urllib.request.urlopen(req, timeout=90) as r: return json.loads(r.read())
            except urllib.error.HTTPError as e:
                if e.code == 401 and attempt == 0: self.tok = None; continue
                sys.exit(f"API {path} 返回 HTTP {e.code}，已停止。")
            except Exception as e:
                print(f"  网络错误（{type(e).__name__}），重试", flush=True); time.sleep(10 * (attempt + 1))
        sys.exit(f"API {path} 多次网络错误，已停止。重跑会续传。")
    def search(self, tile, rel):
        js = self.call("search", {"id": None, "tileId": tile, "releases": [rel], "levels": ["L3"], "productType": "ORTHO-UP"})
        if not isinstance(js, dict): sys.exit("搜索返回格式异常，已停止。")
        if not js.get("status"):
            msg = str(js.get("message", ""))
            if re.search(r"invalid tile", msg, re.I): return None, []
            sys.exit(f"搜索 {tile} {rel} 失败：{msg[:200]}")
        hits = [h for h in js.get("hits", []) if h.get("productType", "ORTHO-UP") == "ORTHO-UP"]
        return js.get("id"), hits

def zip_ok(p):
    try:
        with zipfile.ZipFile(p) as z: return any(n.endswith(".csv") for n in z.namelist())
    except Exception: return False

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def download(sess, tile, rel, fn):
    dst = os.path.join(OUTDIR, fn); part = dst + ".part"
    for attempt in range(4):
        qid, hits = sess.search(tile, rel)  # fresh query id each attempt (ids expire when idle)
        if not any(h["filename"] == fn for h in hits): return "vanished", 0, ""
        u = f"{API}/download/{urllib.parse.quote(fn)}?id={urllib.parse.quote(str(qid), safe='')}"
        try:
            req = urllib.request.Request(u, headers=UA)
            with urllib.request.urlopen(req, timeout=120) as r, open(part, "wb") as f:
                total = int(r.headers.get("Content-Length") or 0); got = 0; nxt = 0.25
                while True:
                    b = r.read(1 << 20)
                    if not b: break
                    f.write(b); got += len(b)
                    if total and got / total >= nxt: print(f"    {fn}: {got/1e6:.0f}/{total/1e6:.0f} MB", flush=True); nxt += 0.25
            if zip_ok(part): os.replace(part, dst); return "ok", os.path.getsize(dst), sha256(dst)
            print(f"    {fn}: 文件不完整，重试", flush=True)
        except urllib.error.HTTPError as e:
            print(f"    {fn}: HTTP {e.code}，重试", flush=True); time.sleep(5)
        except Exception as e:
            print(f"    {fn}: 网络错误（{type(e).__name__}），重试", flush=True); time.sleep(10 * (attempt + 1))
    return "failed", 0, ""

def best(hits):
    return max(hits, key=lambda h: (int(h.get("sub_release") or 0), h["filename"])) if hits else None

def probe(sess, var):
    rels = sess.call("releases")
    lines = [f"{now()} 使用环境变量：{var}（内容不记录）", "访问令牌：获取成功",
             f"存档中的版本（/releases）：{rels}", f"探测瓦片：{PROBE_TILE}"]
    found = {}
    for rel in (rels if isinstance(rels, list) else []):
        if not re.fullmatch(r"\d{4}-\d{4}", str(rel)): continue
        _, hits = sess.search(PROBE_TILE, rel); b = best(hits); found[rel] = b
        lines.append(f"  {rel}: " + (f"{b['filename']}，{int(b.get('filesize') or 0)/1e6:.1f} MB，sub_release={b.get('sub_release')}" if b else "无"))
    open(PROBE_LOG, "w", encoding="utf-8").write("\n".join(lines) + "\n"); print("\n".join(lines))
    for rel in (RA, RB):
        if not found.get(rel): sys.exit(f"{PROBE_TILE} 查不到 {rel}，已停止。请把 probe_log.txt 发给 Claude。")
    print(f"{RA} 与 {RB} 都在，开始下载。", flush=True)

def local_index():
    idx = {}
    for root in (RAW_CACHE, OUTDIR):
        for dp, _, fs in os.walk(root):
            for f in fs:
                if f.startswith("EGMS_L3_") and f.endswith(".zip"): idx.setdefault(f, os.path.join(dp, f))
    return idx

def run_download(sess):
    os.makedirs(OUTDIR, exist_ok=True)
    rows = list(csv.DictReader(open(ORDER, encoding="utf-8"))); idx = local_index()
    new = not os.path.exists(LOG); lf = open(LOG, "a", newline="", encoding="utf-8"); w = csv.writer(lf)
    if new: w.writerow(["utc", "stratum", "rank", "tile", "release", "file", "action", "bytes", "sha256", "path"])
    def get(st, rank, t, rel, hit):
        fn = hit["filename"]
        if fn in idx and zip_ok(idx[fn]):
            w.writerow([now(), st, rank, t, rel, fn, "use_local", os.path.getsize(idx[fn]), "", idx[fn]])
        else:
            res, nb, h = download(sess, t, rel, fn)
            w.writerow([now(), st, rank, t, rel, fn, "downloaded" if res == "ok" else "download_" + res, nb, h, os.path.join(OUTDIR, fn) if res == "ok" else ""])
            if res == "ok": idx[fn] = os.path.join(OUTDIR, fn)
            else: print(f"    {fn}: 未完成（{res}），重跑 run.bat 会续传", flush=True)
        lf.flush()
    for st in ["Al", "S", "C", "Nt"]:
        have = 0
        for r in sorted([r for r in rows if r["stratum"] == st], key=lambda r: int(r["rank"])):
            if have >= QUOTA[st] or int(r["rank"]) > MAX_RANK: break
            t = r["tile"]; plan = {}
            for rel in (RA, RB):
                _, hits = sess.search(t, rel); plan[rel] = best(hits)
            miss = [rel for rel in (RA, RB) if plan[rel] is None]
            if miss:
                w.writerow([now(), st, r["rank"], t, "/".join(miss), "", "skip_not_in_archive", "", "", ""]); lf.flush()
                print(f"[{st} #{r['rank']}] {t}: 存档中缺 {'/'.join(miss)}，跳过", flush=True); continue
            print(f"[{st} #{r['rank']}] {t}: 入选（本层第 {have+1}/{QUOTA[st]} 个）", flush=True)
            for rel in (RA, RB): get(st, r["rank"], t, rel, plan[rel])
            have += 1
        print(f"== {st} 层：{have}/{QUOTA[st]}", flush=True)
    otiles = list(csv.DictReader(open(MANIFEST_O, encoding="utf-8")))
    for k, r in enumerate(otiles, 1):
        t = r["tile"]
        for rel in (RB, RA):
            _, hits = sess.search(t, rel); b = best(hits)
            if b is None:
                w.writerow([now(), "O'", k, t, rel, "", "missing_in_archive", "", "", ""]); lf.flush()
                print(f"[O' {k}/136] {t}: 存档中缺 {rel}", flush=True); continue
            if rel == RA and os.path.basename(r["pathB"]) == b["filename"]:
                w.writerow([now(), "O'", k, t, rel, b["filename"], "use_local_P1B", "", "", r["pathB"]]); lf.flush(); continue
            print(f"[O' {k}/136] {t} {rel}", flush=True); get("O'", k, t, rel, b)
    lf.close(); print("全部结束。日志：download_log.csv")

def fetch_one(sess, tile, rel):
    """Download a single tile/release into tools/egms_download/check/ (structural checks only)."""
    global OUTDIR
    OUTDIR = os.path.join(HERE, "check"); os.makedirs(OUTDIR, exist_ok=True)
    _, hits = sess.search(tile, rel); b = best(hits)
    if not b: sys.exit(f"{tile} {rel}: 存档中没有")
    res, nb, h = download(sess, tile, rel, b["filename"])
    print(f"{b['filename']}: {res}, {nb/1e6:.1f} MB, sha256={h}")

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    var, key = load_key(); sess = Session(key)
    if mode == "fetch": fetch_one(sess, sys.argv[2], sys.argv[3]); sys.exit(0)
    if mode in ("probe", "all"): probe(sess, var)
    if mode in ("download", "all"): run_download(sess)
