#!/usr/bin/env python3
"""Step 5: verify the expected-unchanged list of CORRECTION_PLAN_P84V2.md (READ-ONLY).
  * frozen files: sha256 of protocols.json and complex_metadata.tsv equal the values in the plan
  * Figure 2 panel A inputs (TRAV/TRBV counts) - complex_metadata, and the p84v2 meta agrees with it
  * Figure 2 panel E inputs and duplicate-family leakage per fold (both definitions) - recomputed and
    compared with results/tables/leakage_audit_dual.tsv
  * tcr70 maximum train/test identity - recomputed from the frozen tcr_seq and compared with
    results/tables/percomplex_maxsim.tsv; also, for information, on TCR strings rebuilt from the
    p84v2 rows (every residue, insertions included), which the frozen splits do not use
  * md/: no file under md/ differs from master on this branch
Writes results/review/rebuild_p84v2/step5_verification.tsv.
usage: python scripts/53_verify_unchanged_p84v2.py
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Levenshtein import ratio

sys.path.insert(0, ".")
from src.data import p84v2 as P  # noqa: E402

PLAN = {"data/splits_v2/protocols.json": "3e48a094e94d9758c54187ea924e23caf6c1980bd2d5eb812e85546deb084dca",
        "data/processed/complex_metadata.tsv": "c14641f75b25dc264b2655791e66d2e59492b0468acc606c091fce65c3b4c01f"}
rows = []
def add(item, expected, observed, ok): rows.append(dict(item=item, expected=expected, observed=observed, unchanged=bool(ok)))

for f, h in PLAN.items():
    got = hashlib.sha256(Path(f).read_bytes()).hexdigest()
    add(f"sha256 {f}", h[:16] + "…", got[:16] + "…", got == h)

md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t"); PDBS = list(md.pdb)
proto = json.load(open("data/splits_v2/protocols.json"))
m84 = pd.read_parquet("data/processed/pair_meta_p84v2.parquet").set_index("pdb").loc[PDBS]
ta, tb = md.trav.value_counts(), md.trbv.value_counts()
add("Fig 2A: TRAV/TRBV top-8 counts (complex_metadata)", "61,26,9,7,7,7,6,5 / 49,21,15,15,13,10,8,7",
    ",".join(map(str, ta.head(8).values)) + " / " + ",".join(map(str, tb.head(8).values)),
    list(ta.head(8).values) == [61, 26, 9, 7, 7, 7, 6, 5] and list(tb.head(8).values) == [49, 21, 15, 15, 13, 10, 8, 7])
eq = lambda a, b: np.array([x == y or (pd.isna(x) and pd.isna(y)) for x, y in zip(a, b)])   # one complex lacks TRAV/TRBV in both
ea, eb = eq(m84.trav.values, md.trav.values), eq(m84.trbv.values, md.trbv.values)
add("Fig 2A: p84v2 meta TRAV/TRBV equal complex_metadata (missing == missing)", "174/174", f"{int(ea.sum())}/{int(eb.sum())}", ea.all() and eb.all())

# identities from the frozen tcr_seq
n = len(PDBS); seqs = list(md.tcr_seq)
A = np.eye(n)
for i in range(n):
    for j in range(i + 1, n):
        A[i, j] = A[j, i] = ratio(seqs[i], seqs[j])
idx = {p: i for i, p in enumerate(PDBS)}
pm = pd.read_csv("results/tables/percomplex_maxsim.tsv", sep="\t")
F = []
for name, folds in proto.items():
    for fi, te in enumerate(folds):
        te = [p for p in te if p in idx]; tr = [idx[p] for p in PDBS if p not in set(te)]
        for p in te:
            F.append(dict(protocol=name, fold=fi, pdb=p, ms=A[idx[p], tr].max()))
F = pd.DataFrame(F).merge(pm, on=["protocol", "fold", "pdb"])
t70 = F[F.protocol == "tcr70"]
add("tcr70 maximum train/test identity (frozen tcr_seq)", "0.696", f"{t70.ms.max():.3f}", round(t70.ms.max(), 3) == 0.696)
add("tcr70 test complexes > 0.90", "0 %", f"{100 * (t70.ms > 0.9).mean():.1f} %", (t70.ms > 0.9).sum() == 0)
add("max-identity rows equal percomplex_maxsim.tsv (504)", "504, max |diff| 0", f"{len(F)}, max |diff| {np.abs(F.ms - F.max_sim).max():.1e}",
    len(F) == 504 and np.abs(F.ms - F.max_sim).max() < 1e-12)

# duplicate-family leakage, both definitions (dropna=False reproduces the stored families)
cp = md.cdr3a + "|" + md.cdr3b
fam_s = pd.DataFrame(dict(c=cp, a=md.allele, p=md.peptide)).groupby(["c", "a", "p"], dropna=False).ngroup().values
fam_c = pd.DataFrame(dict(c=cp, p=md.peptide)).groupby(["c", "p"], dropna=False).ngroup().values
au = pd.read_csv("results/tables/leakage_audit_dual.tsv", sep="\t")
mism = 0; t70leak = []
for r in au.itertuples():
    te = set(proto[r.protocol][r.fold]); mt = np.array([p in te for p in PDBS])
    ls, lc = len(set(fam_s[mt]) & set(fam_s[~mt])), len(set(fam_c[mt]) & set(fam_c[~mt]))
    mism += (ls != r.leak_stored) + (lc != r.leak_conservative)
    if r.protocol == "tcr70": t70leak.append(lc)
add("duplicate-family leakage per fold, both definitions, vs leakage_audit_dual.tsv (19 folds)", "0 mismatches", f"{mism} mismatches", mism == 0)
add("tcr70 split duplicate families per fold (conservative)", "0,0,0,0,0", ",".join(map(str, t70leak)), t70leak == [0] * 5)
au38 = au[~((au.protocol == "allele") & (au.n_test < 3))]
e = au38.groupby("protocol").leak_conservative.mean().round(2).to_dict()
add("Fig 2E displayed means (allele fold of n = 2 dropped, as 38b)", "random 7.4, double 0.67, allele 1.25, temporal 0, tcr70 0",
    ", ".join(f"{k} {v:g}" for k, v in e.items()), e == {"allele": 1.25, "double": 0.67, "random": 7.4, "tcr70": 0.0, "temporal": 0.0})

# information only: identity on strings rebuilt from the p84v2 rows (all residues, insertions included)
z = np.load(P.NPZ["4p5A"], allow_pickle=True)
s84 = ["".join(z["tcr_aa"][i]) for i in range(n)]
order = [list(z["pdb"]).index(p) for p in PDBS]; s84 = [s84[i] for i in order]
B = np.eye(n)
for i in range(n):
    for j in range(i + 1, n):
        B[i, j] = B[j, i] = ratio(s84[i], s84[j])
ms84 = [B[idx[p], [idx[q] for q in PDBS if q not in set(fold)]].max() for fold in proto["tcr70"] for p in fold if p in idx]
add("[information] tcr70 max train/test identity on p84v2-row strings (insertions included; not what the splits use)",
    "n/a", f"{max(ms84):.3f}; > 0.70: {sum(x > 0.70 for x in ms84)}; strings changed vs tcr_seq: {sum(a != b for a, b in zip(s84, seqs))}", True)
reg = [str(r) for r in z["row_imgt_region"]]
s84x = ["".join(a for a, rg in zip(z["tcr_aa"][i], reg) if a and rg != "FR4") for i in order]
Bx = np.eye(n)
for i in range(n):
    for j in range(i + 1, n):
        Bx[i, j] = Bx[j, i] = ratio(s84x[i], s84x[j])
msx = [Bx[idx[p], [idx[q] for q in PDBS if q not in set(fold)]].max() for fold in proto["tcr70"] for p in fold if p in idx]
add("[information] same, p84v2-row strings with FR4 118-124 excluded", "n/a", f"{max(msx):.3f}; > 0.70: {sum(x > 0.70 for x in msx)}", True)

diff = subprocess.run(["git", "diff", "--stat", "master", "--", "md"], capture_output=True, text=True).stdout.strip()
add("md/: files differing from master on this branch", "none", diff or "none", diff == "")
T = pd.DataFrame(rows); T.to_csv("results/review/rebuild_p84v2/step5_verification.tsv", sep="\t", index=False)
with pd.option_context("display.width", 250, "display.max_colwidth", 110):
    print(T.to_string(index=False))
print("ALL UNCHANGED" if T.unchanged.all() else "SOMETHING MOVED - see table")
