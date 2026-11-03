#!/usr/bin/env python3
"""Identity bound, Parts A-C (READ-ONLY; data/splits_v2/protocols.json is only read).

Identity = Levenshtein.ratio, imported from scripts/22_build_metadata.py itself.
Four TCR string definitions, all over the six slot segments in SEG order (CDR1/2 IMGT 27-38, 56-65;
CDR3 slot block 105-124):
  (1) frozen      complex_metadata.tcr_seq: one residue per IMGT number (last in file), FR4 included
  (2) full+FR4    every residue (insertion residues restored, IMGT order), FR4 118-124 included
  (3) full-FR4    every residue, FR4 excluded (CDR1/2/3-IMGT only)
  (4) frozen-FR4  definition (1) without FR4
(2) and (3) come from the p84v2 per-row residues (data/processed/pair_maps_p84v2.npz, sha256-checked
against results/review/rebuild_p84v2/artifacts_sha256.tsv). (4) is rebuilt from the same residues with
the old collapse rule; with FR4 kept that rebuild must equal the stored tcr_seq for all 174.
Outputs (results/review/identity_bound/): maxsim_four_definitions.tsv, over070_def2.tsv,
decomposition.tsv, bound_summary.tsv, strings.tsv
usage: python scripts/54a_identity_bound.py
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from src.data import p84v2 as P  # noqa: E402

spec = importlib.util.spec_from_file_location("m22", "scripts/22_build_metadata.py")
m22 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m22)      # main() not called
ratio = m22.ratio
OUT = Path("results/review/identity_bound"); OUT.mkdir(parents=True, exist_ok=True)

npz = P.NPZ["4p5A"]
rec = pd.read_csv("results/review/rebuild_p84v2/artifacts_sha256.tsv", sep="\t").set_index("file")
assert hashlib.sha256(npz.read_bytes()).hexdigest() == rec.loc["data/processed/pair_maps_p84v2.npz", "sha256"], "p84v2 npz hash"
md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t"); PDBS = list(md.pdb); n = len(PDBS)
clu = dict(zip(md.pdb, md.clu70))
proto = json.load(open("data/splits_v2/protocols.json"))
z = np.load(npz, allow_pickle=True)
zi = {str(p): i for i, p in enumerate(z["pdb"])}
seg = [str(s) for s in z["row_segment"]]; num = z["row_imgt_number"].astype(int)
ic = [str(s) for s in z["row_insertion_code"]]; reg = [str(s) for s in z["row_imgt_region"]]
aa = z["tcr_aa"]


def full(p, fr4):
    r = aa[zi[p]]
    return "".join(a for k, a in enumerate(r) if a and (fr4 or reg[k] != "FR4"))


def frozen(p, fr4):
    r = aa[zi[p]]; keep = {}
    for k, a in enumerate(r):
        if not a or (not fr4 and reg[k] == "FR4"): continue
        key = (seg[k], num[k]); rank = P.imgt_rank(num[k], ic[k])
        if key not in keep or rank > keep[key][0]:          # last in file order = highest IMGT rank
            keep[key] = (rank, k)
    return "".join(r[k] for _, k in sorted(keep.values(), key=lambda t: t[1]))


S = {1: [frozen(p, True) for p in PDBS], 2: [full(p, True) for p in PDBS],
     3: [full(p, False) for p in PDBS], 4: [frozen(p, False) for p in PDBS]}
assert S[1] == list(md.tcr_seq), "rebuilt frozen strings differ from complex_metadata.tcr_seq"
DEF = {1: "frozen (collapsed, FR4 incl.)", 2: "full-residue, FR4 incl.", 3: "full-residue, FR4 excl.", 4: "frozen, FR4 excl."}


def is_subseq(a, b):
    it = iter(b); return all(ch in it for ch in a)


n_ins = {p: int(sum(1 for k, a in enumerate(aa[zi[p]]) if a and ic[k])) for p in PDBS}
collapsed = {p: int(sum(1 for a in aa[zi[p]] if a) - len({(seg[k], num[k]) for k, a in enumerate(aa[zi[p]]) if a})) for p in PDBS}
sub_ok = all(is_subseq(S[1][i], S[2][i]) and len(S[2][i]) - len(S[1][i]) == collapsed[p] for i, p in enumerate(PDBS))
same_no_ins = sum(S[1][i] == S[2][i] for i, p in enumerate(PDBS) if n_ins[p] == 0)
pd.DataFrame({"pdb": PDBS, **{f"def{d}": S[d] for d in S}}).to_csv(OUT / "strings.tsv", sep="\t", index=False)

A = {}
for d in S:
    M = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            M[i, j] = M[j, i] = ratio(S[d][i], S[d][j])
    A[d] = M
stored = np.load("results/tables/tcr_identity_matrix.npy")
idx = {p: i for i, p in enumerate(PDBS)}


def fold_max(protocol):
    rows = []
    for f, fold in enumerate(proto[protocol]):
        te = [p for p in fold if p in idx]; tr = [idx[q] for q in PDBS if q not in set(te)]
        for p in te:
            row = dict(protocol=protocol, pdb=p, fold=f, clu70=int(clu[p]))
            for d in S:
                v = A[d][idx[p], tr]; m = v.max(); hits = [PDBS[tr[k]] for k in np.where(np.isclose(v, m))[0]]
                row[f"max_def{d}"] = m; row[f"nn_def{d}"] = hits[0]; row[f"nn_ties_def{d}"] = len(hits)
            rows.append(row)
    return pd.DataFrame(rows)


T = fold_max("tcr70"); R = fold_max("random")
T.to_csv(OUT / "maxsim_four_definitions.tsv", sep="\t", index=False)

# ---------------- the 12 over 0.70 under (2) ----------------
O = T[T.max_def2 > 0.70].sort_values("max_def2", ascending=False).copy()
O["nn_clu70"] = [int(clu[q]) for q in O.nn_def2]
O["pair"] = ["within one cluster" if a == b else "across two clusters" for a, b in zip(O.clu70, O.nn_clu70)]
for d in (1, 3, 4):
    O[f"def{d}_to_def2_neighbour"] = [A[d][idx[p], idx[q]] for p, q in zip(O.pdb, O.nn_def2)]
O.to_csv(OUT / "over070_def2.tsv", sep="\t", index=False)

# ---------------- Part B: decomposition (1) -> (2) for the 12 ----------------
dec = []
for r in O.itertuples():
    t, q1, q2 = idx[r.pdb], idx[r.nn_def1], idx[r.nn_def2]
    id1, id2 = A[1][t, q1], A[2][t, q2]
    sw = A[1][t, q2] - id1                                         # neighbour switch, measured on frozen strings
    mid_t = ratio(S[2][t], S[1][q2])                              # restore test-side insertions only
    mid_n = ratio(S[1][t], S[2][q2])                              # restore neighbour-side insertions only
    dec.append(dict(pdb=r.pdb, nn_def1=r.nn_def1, nn_def2=r.nn_def2, id_def1=id1, id_def2=id2, delta=id2 - id1,
                    neighbour_switch=sw,
                    restore_test_then_nn=f"{mid_t - A[1][t, q2]:+.4f} then {id2 - mid_t:+.4f}",
                    restore_nn_then_test=f"{mid_n - A[1][t, q2]:+.4f} then {id2 - mid_n:+.4f}",
                    restoration_total=id2 - A[1][t, q2],
                    restored_residues_test=len(S[2][t]) - len(S[1][t]), restored_residues_nn=len(S[2][q2]) - len(S[1][q2]),
                    fr4_effect_same_neighbour=A[2][t, q2] - A[3][t, q2], id_def3=T.loc[T.pdb == r.pdb, "max_def3"].iloc[0],
                    over070_def1=id1 > 0.70, over070_def3=bool(T.loc[T.pdb == r.pdb, "max_def3"].iloc[0] > 0.70),
                    over070_def4=bool(T.loc[T.pdb == r.pdb, "max_def4"].iloc[0] > 0.70)))
DEC = pd.DataFrame(dec); DEC.to_csv(OUT / "decomposition.tsv", sep="\t", index=False)

# ---------------- Part C: the bound under each definition ----------------
summ = []
for name, F in (("tcr70", T), ("random", R)):
    for d in S:
        v = F[f"max_def{d}"].values
        summ.append(dict(protocol=name, definition=f"({d}) {DEF[d]}", n=len(v), max=v.max(), p95=np.quantile(v, 0.95),
                         median=np.median(v), mean=v.mean(), n_gt_070=int((v > 0.70).sum()), n_gt_080=int((v > 0.80).sum()),
                         n_gt_090=int((v > 0.90).sum()), n_ge_090=int((v >= 0.90).sum())))
SUM = pd.DataFrame(summ); SUM.to_csv(OUT / "bound_summary.tsv", sep="\t", index=False)

pm = pd.read_csv("results/tables/percomplex_maxsim.tsv", sep="\t")
chk = T.merge(pm, on=["pdb", "fold"]).query("protocol_x == 'tcr70' and protocol_y == 'tcr70'")
with pd.option_context("display.width", 250, "display.max_colwidth", 60, "display.float_format", "{:.4f}".format):
    print(f"checks: def(1) rebuilt == stored tcr_seq 174/174; def(1) matrix == tcr_identity_matrix.npy (max |diff| {np.abs(A[1] - stored).max():.1e}); "
          f"def(1) tcr70 maxima == percomplex_maxsim (max |diff| {np.abs(chk.max_def1 - chk.max_sim).max():.1e}); "
          f"def(1) is a subsequence of def(2) with exactly the restored residues added: {sub_ok}; "
          f"complexes without insertions: def(1) == def(2) in {same_no_ins}/{sum(1 for p in PDBS if n_ins[p] == 0)}")
    print("\nPart A - over 0.70 under (2):"); print(O[["pdb", "fold", "clu70", "nn_def2", "nn_clu70", "pair", "max_def1", "max_def2", "max_def3", "max_def4",
                                                         "def1_to_def2_neighbour", "def3_to_def2_neighbour", "def4_to_def2_neighbour", "nn_def1", "nn_def3", "nn_def4"]].to_string(index=False))
    print("\nPart B - decomposition:"); print(DEC.to_string(index=False))
    print("\nPart C - bound:"); print(SUM.to_string(index=False))
