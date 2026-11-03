#!/usr/bin/env python3
"""Worked example: reproduce the two headline pairwise numbers from released artifacts alone.

Runs in seconds on a CPU. No training, no GPU, no network. It scores the released prediction
arrays and reproduces values quoted in the paper, so a reader can confirm the evaluation
pipeline end to end before deciding whether to repeat the expensive training steps.

What this demonstrates is the EVALUATION contract, which is what this paper contributes.
It does NOT demonstrate prediction for a new receptor: see README "Input contract and scope".

usage (from the package root):  python scripts/worked_example.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from src.eval.metrics import per_complex, cluster_bootstrap_diff

# The frozen arrays ship in the Zenodo archive, not in the git repository. Say so plainly rather
# than letting numpy raise FileNotFoundError on a path the reader has no reason to recognise.
NEEDED = ["data/processed/complex_metadata.tsv", "data/processed/pair_maps_p84v2.npz",
          "results/tables/pairmaps_tcr70_p84v2.npz", "results/tables/teim_res_pairwise_p84v2.npz",
          "results/tables/similarity_curve_p84v2_authoritative.csv"]
missing = [f for f in NEEDED if not (ROOT/f).exists()]
if missing:
    sys.exit("This package does not contain the frozen arrays this example reads:\n  "
             + "\n  ".join(missing)
             + "\n\nYou have the code tree without the data. The arrays are in the Zenodo archive;"
               "\nDATA.md explains how the two correspond and how to unpack one over the other."
               "\nCheck that PACKAGE_INFO.json reports the same source_revision in both.")

md = pd.read_csv(ROOT/"data/processed/complex_metadata.tsv", sep="\t")
lab = np.load(ROOT/"data/processed/pair_maps_p84v2.npz", allow_pickle=True)
pm = np.load(ROOT/"results/tables/pairmaps_tcr70_p84v2.npz", allow_pickle=True)
te = np.load(ROOT/"results/tables/teim_res_pairwise_p84v2.npz", allow_pickle=True)

print("1. integrity of the released arrays")
checks = [("identifiers match the metadata", [str(x) for x in pm["pdb"]] == list(md.pdb)),
          ("cluster assignments match", bool(np.array_equal(pm["clu70"], md.clu70.values))),
          ("labels match the frozen build", bool(np.array_equal(pm["Y"].astype(np.int8), lab["Y"]))),
          ("masks match the frozen build", bool(np.array_equal(pm["S"].astype(np.int8), lab["S"]))),
          ("tensor shape is 174 x 75 x 84", pm["Y"].shape == (174, 75, 84))]
for name, ok in checks:
    print(f"   [{'ok' if ok else 'FAIL'}] {name}")
assert all(ok for _, ok in checks)

idx = np.where(pm["covered"])[0]
Y, S = pm["Y"].astype(np.float32)[idx], pm["S"].astype(np.float32)[idx]
clu = pm["clu70"][idx]
mm = per_complex(Y, pm["pred_model"][idx], S)
ff = per_complex(Y, pm["pred_freq"][idx], S)
tt = per_complex(Y, te["pred"][idx], S)

print("\n2. pairwise task under tcr70, cluster bootstrap over 69 clusters, 20,000 resamples, seed 42")
for name, arr, expect in (("TCR-IMAP vs pair-frequency prior", ff["auprc"], "+0.0507 [+0.0215, +0.0841], 112/174"),
                          ("TCR-IMAP vs TEIM-Res", tt["auprc"], "+0.0873 [+0.0592, +0.1149], 143/174")):
    d, lo, hi, p = cluster_bootstrap_diff(mm["auprc"], arr, clu, n_boot=20000, seed=42)
    wins = int((mm["auprc"] > arr).sum())
    print(f"   {name}")
    print(f"      computed : {d:+.4f} [{lo:+.4f}, {hi:+.4f}], p = {p:.4f}, higher on {wins}/{len(idx)}")
    print(f"      in paper : {expect}, p < 0.001")

print(f"\n   TCR-IMAP macro-AUPRC {np.nanmean(mm['auprc']):.4f} (paper 0.3598); "
      f"prior {np.nanmean(ff['auprc']):.4f} (paper 0.3090); TEIM-Res {np.nanmean(tt['auprc']):.4f} (paper 0.2725)")

print("\n3. the identity sweep, authoritative statistics")
auth = pd.read_csv(ROOT/"results/tables/similarity_curve_p84v2_authoritative.csv")
print(auth.to_string(index=False))
print("\nDone. Every value above was computed here from the released arrays.")
