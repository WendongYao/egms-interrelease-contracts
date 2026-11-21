"""Check every SHA-256 entry in the frozen hash files of the four contracts.

Hash files were written from different working directories, so each entry is
looked up in order: next to the hash file, contract root, code/, results/.
Documented exceptions resolved here (see README.md, section "Hash check"):
  * Confirm01 CONTRACT_CN.md: the frozen hash covers lines 1-121, i.e. the text
    before deviation entries were appended (rule: deviations appended at end).
  * Result files listed as  *_N'.csv / *_O'.csv  are stored as *_Nprime.csv /
    *_Oprime.csv (same bytes).
  * Confirm03 lists inputs from other contracts: data/* from Confirm01 data/,
    /home/claude/c02/code/rsc/* = Confirm02 code/rsc/*.
  * Later code versions recorded in RUN_LOG / deviations are reported as
    SUPERSEDED when the current file matches the later hash given there.
Usage: python verify_hashes.py [contracts_dir]"""
import hashlib, os, re, sys, glob
root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "contracts")
pat = re.compile(r"^([0-9a-f]{64})\s+\*?(.+?)\s*$")
LATER = {  # file -> later hash documented in RUN_LOG_CN.md / CONTRACT_CN.md deviations
    ("Confirm01_JAG", "code/s89_part2.py"): "d6c735c1193cb2704aaa6aad20f437cb59da4ac15b1ebdc61594d71053c964fa",
    ("Confirm01_JAG", "code/s91_manifest_P2.py"): "add1a6c15926828703216d8fc8951e10fec66cabd29daaa3773b40f4554f162f",
}
def sha_bytes(b): return hashlib.sha256(b).hexdigest()
def sha(p):
    with open(p, "rb") as f: return sha_bytes(f.read())
contracts = sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d)))
tot = {}
for c in contracts:
    cd = os.path.join(root, c)
    hfiles = [os.path.join(cd, "README_HASH.txt")] + sorted(glob.glob(os.path.join(cd, "code", "*HASH*.txt"))) + sorted(glob.glob(os.path.join(cd, "results", "*HASH*.txt")))
    for hf in hfiles:
        if not os.path.isfile(hf): continue
        for line in open(hf, encoding="utf-8", errors="replace"):
            m = pat.match(line.strip())
            if not m: continue
            h, p = m.group(1), m.group(2).strip('"')
            q = p.replace("_N'.", "_Nprime.").replace("_O'.", "_Oprime.")
            cands = [os.path.join(b, x) for x in dict.fromkeys([p, q]) for b in (os.path.dirname(hf), cd, os.path.join(cd, "code"), os.path.join(cd, "results"))]
            if p.startswith("/home/claude/c02/code/"):
                cands.append(os.path.join(root, "Confirm02_Mechanism", "code", p.split("/home/claude/c02/code/")[1]))
            if q.startswith("data/"):
                cands += [os.path.join(root, o, q) for o in contracts]
            f = next((x for x in cands if os.path.isfile(x)), None)
            if f is None: st = "ABSENT"
            elif sha(f) == h: st = "OK"
            elif c == "Confirm01_JAG" and p == "CONTRACT_CN.md" and sha_bytes(b"".join(open(f, "rb").readlines()[:121])) == h: st = "OK_PREFIX"
            elif LATER.get((c, p)) == sha(f): st = "SUPERSEDED"
            else: st = "MISMATCH"
            tot[st] = tot.get(st, 0) + 1
            if st != "OK": print(f"{st:10s} {c}/{os.path.relpath(hf, cd)} :: {p}")
print("TOTAL", tot)
