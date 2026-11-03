#!/usr/bin/env python3
"""Composition of the identity-threshold sweep (Figure 3A/B; SI Figure S9): per threshold, the test and training
complexes each fold actually used, and the train/test similarity actually achieved.

Authorised by the author on 26 October 2026 (results/review/SI_FIGURES_REPORT.md, S9). Arithmetic only, on committed
inputs: no training, inference or scoring.

Inputs and the sweep's own rules (scripts/25b_similarity_curve_p82v1.py, executed for p84v2 by
scripts/53_similarity_curve_p84v2.py):
- similarity: results/tables/tcr_identity_matrix.npy (committed; Levenshtein ratio of complex_metadata.tcr_seq, the
  strings the sweep compares), checked here against python-Levenshtein on those strings
- folds: np.random.default_rng(0).permutation(174) split five ways (np.array_split), in pair-map order
- training complexes kept at threshold thr: every non-test complex t with max_j Sim[t, test_j] < thr; a fold is
  used only if at least 30 remain

Per threshold it writes, to results/tables/similarity_curve_composition_p84v2.csv:
- n_test (complexes evaluated) and the per-fold test and training-set sizes;
- the ACHIEVED train/test similarity: for each test complex, its maximum similarity to the training complexes its
  fold actually kept; then the mean and the maximum of that over the test complexes;
- mean_max_sim_column: the existing similarity_curve_p84v2.csv column, reproduced here with the sweep's own formula.
  That column is, for each test complex, the maximum similarity to ANY other complex of the 174 whose similarity is
  below the threshold - including complexes in the same test fold and training complexes the filter removed - so it
  is an upper bound on the achieved value, not the achieved value.
Per-fold rows go to results/tables/similarity_curve_composition_folds_p84v2.csv.
usage: python scripts/57_sweep_composition_p84v2.py
"""
import numpy as np, pandas as pd
from Levenshtein import ratio

T = "results/tables/"
THRESHOLDS = [1.01, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70]          # the sweep's own list
RD = dict(float_precision="round_trip")

md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t")
Z = np.load("data/processed/pair_maps_p84v2.npz", allow_pickle=True)
assert [str(p) for p in Z["pdb"]] == md.pdb.tolist()               # the sweep indexes complexes in this order
Sim = np.load(T + "tcr_identity_matrix.npy")
n = len(md)
assert Sim.shape == (n, n) and np.allclose(np.diag(Sim), 1.0)
seqs = md.tcr_seq.tolist()
R = np.eye(n)
for i in range(n):
    for j in range(i + 1, n):
        R[i, j] = R[j, i] = ratio(seqs[i], seqs[j])
assert np.array_equal(R, Sim), np.abs(R - Sim).max()
print(f"tcr_identity_matrix.npy equals Levenshtein.ratio of complex_metadata.tcr_seq exactly ({n} x {n})")

rng = np.random.default_rng(0)
folds = np.array_split(rng.permutation(n), 5)
assert sorted(np.concatenate(folds).tolist()) == list(range(n)) and [len(f) for f in folds] == [35, 35, 35, 35, 34]

sc = pd.read_csv(T + "similarity_curve_p84v2.csv", **RD).set_index("threshold")
rows, frows = [], []
for thr in THRESHOLDS:
    used = np.zeros(n, bool)
    ach = np.full(n, np.nan)
    for k, te in enumerate(folds):
        tr = np.setdiff1d(np.arange(n), te)
        keep = np.array([Sim[t, te].max() < thr for t in tr])        # the sweep's rule, verbatim
        tr2 = tr[keep]
        frows.append(dict(threshold=thr, fold=k, n_test=len(te), n_train_available=len(tr), n_train_kept=len(tr2),
                          n_train_removed=len(tr) - len(tr2), fold_used=len(tr2) >= 30))
        if len(tr2) < 30:
            continue
        assert not np.isin(tr2, te).any()
        used[te] = True
        ach[te] = Sim[np.ix_(te, tr2)].max(axis=1)
        assert (ach[te] < thr).all()
    idx = np.where(used)[0]
    # the existing column, with the sweep's own expression
    ub = np.mean([Sim[i, np.setdiff1d(np.arange(n), [i])][Sim[i, np.setdiff1d(np.arange(n), [i])] < thr].max()
                  if (Sim[i, np.setdiff1d(np.arange(n), [i])] < thr).any() else 0 for i in idx])
    rec = sc.loc[round(thr, 2)]
    assert len(idx) == int(rec.n), (thr, len(idx), rec.n)
    assert round(ub, 4) == rec.mean_max_sim, (thr, ub, rec.mean_max_sim)
    assert (ach[idx] <= np.array([Sim[i, np.setdiff1d(np.arange(n), [i])][Sim[i, np.setdiff1d(np.arange(n), [i])] < thr].max()
                                  for i in idx]) + 1e-15).all()     # achieved never exceeds the column's per-complex value
    fr = [r for r in frows if r["threshold"] == thr]
    rows.append(dict(threshold=thr, n_test=len(idx), folds_used=sum(r["fold_used"] for r in fr),
                     n_train_min=min(r["n_train_kept"] for r in fr), n_train_max=max(r["n_train_kept"] for r in fr),
                     n_train_mean=float(np.mean([r["n_train_kept"] for r in fr])),
                     n_train_available_min=min(r["n_train_available"] for r in fr),
                     n_train_available_max=max(r["n_train_available"] for r in fr),
                     achieved_max_sim_mean=float(np.mean(ach[idx])), achieved_max_sim_max=float(np.max(ach[idx])),
                     achieved_max_sim_median=float(np.median(ach[idx])),
                     mean_max_sim_column=float(ub)))
D = pd.DataFrame(rows)
D.to_csv(T + "similarity_curve_composition_p84v2.csv", index=False)
pd.DataFrame(frows).to_csv(T + "similarity_curve_composition_folds_p84v2.csv", index=False)
print("n test = 174 at every threshold (= similarity_curve_p84v2.csv n); mean_max_sim reproduced at all 7 thresholds")
print(D.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
