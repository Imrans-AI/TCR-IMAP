#!/usr/bin/env python3
"""p84v2 pairwise contact labels: 18b's definition, with residues keyed by (IMGT number, insertion
code) throughout and every residue in its own row/column (results/review/CORRECTION_PLAN_P84V2.md).

Same as 18b: STCRDab IMGT structures (data/structures/stcrdab/imgt), chains from
data/interim/stcrdab_index.tsv, standard residues only (id[0] == " "), heavy atoms (element != "H"),
contact = minimum heavy-atom distance <= CUTOFF (4.5 A; 6.0 A for the variant), same 174 complexes
in the same order (taken from pair_maps_p82v1.npz, read only).
Different from 18b:
  * TCR rows: every (segment, IMGT number, insertion code) occupied in >= 1 of the 174 complexes,
    within the 18b slot ranges (IMGT 27-38, 56-65, 105-124), in IMGT order; never-occupied slots dropped
  * MHC columns: every (IMGT number, insertion code) occupied within IMGT 50-90, 1050-1090
  * no residue is ever replaced by another sharing its number
  * per-row residue identities (tcr_aa, mhc_aa, *_resid) are stored, for position-placed tokens
Outputs (refuses to overwrite): data/processed/pair_maps_p84v2[_6A].npz, pair_meta_p84v2[_6A].parquet
usage: python scripts/53_build_pair_labels_p84v2.py [--cutoff 6.0]
"""
import argparse
import datetime
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from scipy.spatial.distance import cdist

sys.path.insert(0, ".")
from src.data.p84v2 import imgt_rank, SEG_ORDER, SLOT_START, SLOT_LENS_P82V1, POS82  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--cutoff", type=float, default=4.5); a = ap.parse_args()
CUTOFF = a.cutoff
TAG = "" if CUTOFF == 4.5 else f"_{CUTOFF:g}A"
OUT_NPZ = Path(f"data/processed/pair_maps_p84v2{TAG}.npz"); OUT_META = Path(f"data/processed/pair_meta_p84v2{TAG}.parquet")
for f in (OUT_NPZ, OUT_META):
    assert not f.exists(), f"refusing to overwrite existing {f}"

def region(p):
    return "CDR1" if 27 <= p <= 38 else "CDR2" if 56 <= p <= 65 else "CDR3" if 105 <= p <= 117 else "FR4"

pdbs = [str(p) for p in np.load("data/processed/pair_maps_p82v1.npz", allow_pickle=True)["pdb"]]
idx = pd.read_csv("data/interim/stcrdab_index.tsv", sep="\t").drop_duplicates("pdb").set_index("pdb")
parser = PDBParser(QUIET=True)

# ---- pass 1: residues per complex keyed by (number, icode); the axes are the union over complexes
per = {}
rowkeys, colkeys = set(), set()
for pdb in pdbs:
    r = idx.loc[pdb]
    model = parser.get_structure(pdb, f"data/structures/stcrdab/imgt/{pdb}.pdb")[0]
    tcr = {}
    for tag, ch in (("a", r.Achain), ("b", r.Bchain)):
        for x in model[str(ch)]:
            if x.id[0] != " ": continue
            for seg, lo in SLOT_START.items():
                key = f"{seg}{tag}"
                if lo <= x.id[1] < lo + SLOT_LENS_P82V1[key]:
                    k = (key, x.id[1], x.id[2].strip())
                    assert k not in tcr, f"{pdb}: duplicate residue key {k}"
                    tcr[k] = x; rowkeys.add(k)
    mhc = {}
    for x in model[str(r.mhc_chain1)]:
        if x.id[0] == " " and x.id[1] in set(POS82):
            k = (x.id[1], x.id[2].strip())
            assert k not in mhc, f"{pdb}: duplicate MHC key {k}"
            mhc[k] = x; colkeys.add(k)
    per[pdb] = (tcr, mhc, r)

rows = sorted(rowkeys, key=lambda k: (SEG_ORDER.index(k[0]), imgt_rank(k[1], k[2])))
cols = sorted(colkeys, key=lambda k: imgt_rank(*k))
R, C = {k: i for i, k in enumerate(rows)}, {k: j for j, k in enumerate(cols)}

# ---- pass 2: labels
n = len(pdbs)
Y = np.zeros((n, len(rows), len(cols)), np.int8); S = np.zeros_like(Y)
tcr_aa = np.full((n, len(rows)), "", dtype="<U1"); mhc_aa = np.full((n, len(cols)), "", dtype="<U1")
tcr_resid = np.full((n, len(rows)), "", dtype="<U6"); mhc_resid = np.full((n, len(cols)), "", dtype="<U6")
meta = []
for c, pdb in enumerate(pdbs):
    tcr, mhc, r = per[pdb]
    ti = sorted(tcr, key=lambda k: R[k]); mj = sorted(mhc, key=lambda k: C[k])
    for k in ti:
        tcr_aa[c, R[k]] = seq1(tcr[k].get_resname()); tcr_resid[c, R[k]] = f"{k[1]}{k[2]}"
    for k in mj:
        mhc_aa[c, C[k]] = seq1(mhc[k].get_resname()); mhc_resid[c, C[k]] = f"{k[0]}{k[1]}"
    ta, to, ma, mo = [], [], [], []
    for i, k in enumerate(ti):
        for at in tcr[k]:
            if at.element != "H": ta.append(at.coord); to.append(i)
    for j, k in enumerate(mj):
        for at in mhc[k]:
            if at.element != "H": ma.append(at.coord); mo.append(j)
    hit = cdist(np.array(ta, float), np.array(ma, float)) <= CUTOFF
    M = np.zeros((len(ti), len(mj)), bool)
    ii, jj = np.nonzero(hit); M[np.array(to)[ii], np.array(mo)[jj]] = True
    rr = np.array([R[k] for k in ti]); cc = np.array([C[k] for k in mj])
    S[c][np.ix_(rr, cc)] = 1
    Y[c][np.ix_(rr, cc)] = M.astype(np.int8)
    meta.append(dict(pdb=pdb, allele=r.allele if isinstance(r.allele, str) else "NA", trbv=r.beta_subgroup, trav=r.alpha_subgroup))

git = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
np.savez_compressed(
    OUT_NPZ, Y=Y, S=S, pdb=np.array(pdbs),
    row_label=np.array([f"{s}:{p}{i}" for s, p, i in rows]), row_segment=np.array([s for s, _, _ in rows]),
    row_imgt_number=np.array([p for _, p, _ in rows], np.int64), row_insertion_code=np.array([i for _, _, i in rows]),
    row_imgt_region=np.array([region(p) for _, p, _ in rows]),
    col_label=np.array([f"{p}{i}" for p, i in cols]), col_imgt_number=np.array([p for p, _ in cols], np.int64),
    col_insertion_code=np.array([i for _, i in cols]),
    tcr_aa=tcr_aa, mhc_aa=mhc_aa, tcr_resid=tcr_resid, mhc_resid=mhc_resid,
    cutoff=np.array(CUTOFF), candidate_set=np.array(f"imgt84_p84v2{TAG}"),
    provenance=np.array(f"scripts/53_build_pair_labels_p84v2.py | git {git} | built {datetime.datetime.now().isoformat(timespec='seconds')} | "
                        f"structures data/structures/stcrdab/imgt | chains data/interim/stcrdab_index.tsv | complex set+order from "
                        f"data/processed/pair_maps_p82v1.npz | keying (IMGT number, insertion code) | heavy atoms element!=H | "
                        f"contact min distance <= {CUTOFF} A | plan results/review/CORRECTION_PLAN_P84V2.md"))
pd.DataFrame(meta).to_parquet(OUT_META, index=False)
print(f"wrote {OUT_NPZ}: Y/S {Y.shape}; resolved {int(S.sum())}; positives {int(Y.sum())} ({Y.sum() / S.sum():.6f}); cutoff {CUTOFF}")
print(f"wrote {OUT_META}")
