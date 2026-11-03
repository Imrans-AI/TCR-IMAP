#!/usr/bin/env python3
"""Bulk-download STCRDab IMGT-numbered structures + per-entry summaries."""
import argparse, re, sys, time
from pathlib import Path
import requests

BASE = "https://opig.stats.ox.ac.uk/webapps/stcrdab-stcrpred"
LISTING = f"{BASE}/Browser?all=true"
HEADERS = {"User-Agent": "tcr-imap-research/0.1 (academic use)"}

def get_pdb_ids(s):
    r = s.get(LISTING, headers=HEADERS, timeout=120); r.raise_for_status()
    ids = re.findall(r"/stcrdab-stcrpred/pdb/([0-9a-zA-Z]{4})", r.text)
    seen, out = set(), []
    for i in ids:
        i = i.lower()
        if i not in seen:
            seen.add(i); out.append(i)
    return out

def fetch(s, url, dest, sleep):
    if dest.exists() and dest.stat().st_size > 0:
        return "skip"
    try:
        r = s.get(url, headers=HEADERS, timeout=60)
        if r.status_code != 200 or not r.content.strip():
            return f"fail({r.status_code})"
        dest.write_bytes(r.content); time.sleep(sleep); return "ok"
    except Exception as e:
        return f"error({type(e).__name__})"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/structures/stcrdab")
    ap.add_argument("--sleep", type=float, default=0.4)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    out = Path(a.out)
    (out/"imgt").mkdir(parents=True, exist_ok=True)
    (out/"summary").mkdir(parents=True, exist_ok=True)
    s = requests.Session()
    print("Fetching listing ...", flush=True)
    ids = get_pdb_ids(s)
    if not ids: sys.exit("No PDB IDs found - page layout may have changed.")
    print(f"Found {len(ids)} entries.", flush=True)
    if a.limit: ids = ids[:a.limit]
    fails = []
    for n, p in enumerate(ids, 1):
        s1 = fetch(s, f"{BASE}/pdb/{p}", out/"imgt"/f"{p}.pdb", a.sleep)
        s2 = fetch(s, f"{BASE}/summary/{p}", out/"summary"/f"{p}.tsv", a.sleep)
        if "ok" not in (s1,s2) and "skip" not in (s1,s2): fails.append((p,s1,s2))
        if n % 25 == 0 or n == len(ids):
            print(f"  {n}/{len(ids)} last={p} pdb={s1} summary={s2}", flush=True)
    (out/"FETCH_MANIFEST.txt").write_text(
        f"source: {LISTING}\nfetched: {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"entries_listed: {len(ids)}\n"
        f"pdb_files: {len(list((out/'imgt').glob('*.pdb')))}\n"
        f"summary_files: {len(list((out/'summary').glob('*.tsv')))}\n"
        f"failures: {len(fails)}\n" + "".join(f"  {x} {y} {z}\n" for x,y,z in fails))
    print("\nDone. Manifest written.")
    if fails: print(f"{len(fails)} entries had issues - re-run to retry.")

if __name__ == "__main__":
    main()
