#!/usr/bin/env python3
"""TCR residue graphs for Struct-TCS under the p84v2 layout: 20_build_graphs.py's definition
(C-beta, C-alpha for Gly, within 8 A; self-loops; nodes without coordinates isolated), with nodes =
the 75 p84v2 rows keyed by (IMGT number, insertion code) instead of 84 number-only slots.

Validation first (--validate, default on): the same code run in p82v1 mode (number-only keying,
last residue in file wins, 84 slots) must reproduce data/processed/tcr_graphs.npz exactly; only
then is the p84v2 graph written.
Output (refuses to overwrite): data/processed/tcr_graphs_p84v2.npz (A: 174 x 75 x 75, pdb, row_label)
usage: python scripts/53_build_graphs_p84v2.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

sys.path.insert(0, ".")
from src.data.p84v2 import load, SLOT_START, SLOT_LENS_P82V1, SEG_ORDER  # noqa: E402

OUT = Path("data/processed/tcr_graphs_p84v2.npz")
assert not OUT.exists(), f"refusing to overwrite existing {OUT}"
CUT = 8.0
d = load()
PDBS = d["pdb"]
rows = list(zip([lab.split(":")[0] for lab in d["row_label"]], d["row_imgt_number"], d["row_insertion_code"]))
R = {k: i for i, k in enumerate(rows)}
idx = pd.read_csv("data/interim/stcrdab_index.tsv", sep="\t").set_index("pdb")
parser = PDBParser(QUIET=True)


def cb(res):
    if "CB" in res: return res["CB"].coord
    if "CA" in res: return res["CA"].coord
    return None


def adjacency(coords):
    ok = ~np.isnan(coords[:, 0])
    dd = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    adj = (dd <= CUT) & ok[:, None] & ok[None, :]
    np.fill_diagonal(adj, True)
    return adj.astype(np.float32)


A82 = np.zeros((len(PDBS), 84, 84), np.float32); A84 = np.zeros((len(PDBS), len(rows), len(rows)), np.float32)
off82 = {k: sum(SLOT_LENS_P82V1[j] for j in SEG_ORDER[:SEG_ORDER.index(k)]) for k in SEG_ORDER}
for n, pdb in enumerate(PDBS):
    r = idx.loc[pdb]
    if isinstance(r, pd.DataFrame): r = r.iloc[0]
    m = parser.get_structure(pdb, f"data/structures/stcrdab/imgt/{pdb}.pdb")[0]
    c82 = np.full((84, 3), np.nan, np.float32); c84 = np.full((len(rows), 3), np.nan, np.float32)
    for tag, ch in (("a", r.Achain), ("b", r.Bchain)):
        std = [x for x in m[str(ch)] if x.id[0] == " "]
        tc = {x.id[1]: x for x in std}                                  # p82v1 keying (as 20_build_graphs)
        for seg, lo in SLOT_START.items():
            key = f"{seg}{tag}"
            for p in range(lo, lo + SLOT_LENS_P82V1[key]):
                if p in tc:
                    v = cb(tc[p])
                    if v is not None: c82[off82[key] + (p - lo)] = v
            for x in std:                                               # p84v2 keying
                if lo <= x.id[1] < lo + SLOT_LENS_P82V1[key]:
                    v = cb(x)
                    if v is not None: c84[R[(key, x.id[1], x.id[2].strip())]] = v
    A82[n] = adjacency(c82); A84[n] = adjacency(c84)

ref = np.load("data/processed/tcr_graphs.npz")["A"]
same = bool((A82 == ref).all())
print(f"validation: p82v1-mode graphs == data/processed/tcr_graphs.npz: {same} (mismatching entries {int((A82 != ref).sum())})")
if not same:
    print("STOP: the graph code does not reproduce the p82v1 graphs"); sys.exit(1)
np.savez_compressed(OUT, A=A84, pdb=np.array(PDBS), row_label=np.array(d["row_label"]))
deg = A84.sum(-1)
print(f"wrote {OUT}: {A84.shape}; mean degree over present nodes {deg[deg > 0].mean():.2f} "
      f"(p82v1 {ref.sum(-1)[ref.sum(-1) > 0].mean():.2f})")
