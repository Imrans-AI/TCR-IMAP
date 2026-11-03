#!/usr/bin/env python3
"""Comparator subset: compare the two settings on the SAME 22 test complexes.

The reported comparator-subset result contrasts a single random split (22 test complexes) with
tcr70 (88 test complexes) at each cutoff. Those two cells differ in the partitioning protocol AND
in which complexes are scored, so the contrast confounds the protocol with the test set.

This script removes that confound using only STORED per-complex arrays -- nothing is retrained.
It restricts the tcr70 setting to exactly the 22 complexes the single random split evaluated, and
recomputes the paired difference there, at both the 4.5 A and 6 A cutoffs.

The 88-complex tcr70 result is PRESERVED and reported alongside; it is not replaced.

usage: python scripts/65_same22_comparator.py
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, ".")
from src.eval.metrics import cluster_bootstrap_diff

ROOT = Path(__file__).resolve().parent.parent
CACHE45 = ROOT/"results/tables/_partial/p84v2/31b_structcs_45A"
CACHE6 = ROOT/"results/tables/_partial/p84v2/31b_structcs_6A"
PERSEED = ROOT/"results/tables/_partial/p84v2/31b_perseed"
OUT = ROOT/"results/review/same22"; OUT.mkdir(parents=True, exist_ok=True)
N_BOOT, SEED = 20000, 42
COMPARATORS = ["Struct-TCS", "Seq-TCS"]


def load(cut, setting):
    for d in (CACHE45, CACHE6, PERSEED):
        f = d/f"{cut}_{setting}.npz"
        if f.exists():
            return json.loads(str(np.load(f, allow_pickle=True)["__percomplex__"]))
    raise FileNotFoundError(f"{cut}_{setting}")


def main():
    rows = []
    for cut in ("4p5A", "6A"):
        rnd = load(cut, "their_random")
        t70 = load(cut, "our_tcr70")
        pr, pt = list(rnd["pdb"]), list(t70["pdb"])
        shared = [p for p in pr if p in set(pt)]
        print(f"\n=== cutoff {cut} ===")
        print(f"single random split: {len(pr)} test complexes;  tcr70: {len(pt)};  "
              f"random-split complexes also in the tcr70 evaluation: {len(shared)}")
        if len(shared) != len(pr):
            print(f"  NOT a subset -- {len(pr)-len(shared)} random-split complexes are absent from tcr70: "
                  f"{sorted(set(pr)-set(pt))}")
        ti = {p: i for i, p in enumerate(pt)}
        ri = {p: i for i, p in enumerate(pr)}
        for comp in COMPARATORS:
            for label, src, idx, clu in (
                ("single random split, n = %d" % len(pr), rnd, [ri[p] for p in pr],
                 np.array(rnd["clu70"])),
                ("tcr70, full subset, n = %d" % len(pt), t70, [ti[p] for p in pt],
                 np.array(t70["clu70"])),
                ("tcr70, SAME %d complexes" % len(shared), t70, [ti[p] for p in shared],
                 np.array([t70["clu70"][ti[p]] for p in shared]))):
                a = np.array(src["Ours_auprc"], float)[idx]
                b = np.array(src[f"{comp}_auprc"], float)[idx]
                d, lo, hi, p = cluster_bootstrap_diff(a, b, clu, n_boot=N_BOOT, seed=SEED)
                rows.append(dict(cutoff=cut.replace("p5", ".5").replace("A", " A"),
                                 comparator=comp, setting=label, n=len(idx),
                                 n_clusters=len(np.unique(clu)),
                                 ours=round(float(np.nanmean(a)), 4),
                                 comp=round(float(np.nanmean(b)), 4),
                                 diff=round(float(d), 4), lo=round(float(lo), 4),
                                 hi=round(float(hi), 4), p=round(float(p), 4),
                                 ci_excludes_zero=bool(lo > 0 or hi < 0)))
    t = pd.DataFrame(rows)
    t.to_csv(OUT/"same22_comparison.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    print("\n=== paired differences, TCR-IMAP minus comparator (20,000 cluster resamples, seed 42) ===")
    print(t.to_string(index=False))
    print("\nNote: the 'SAME n' rows isolate the protocol by holding the test complexes fixed; the "
          "full-subset tcr70 rows are preserved and remain the reported comparator-subset result.")
    print(f"\nwrote {OUT/'same22_comparison.tsv'}")


if __name__ == "__main__":
    main()
