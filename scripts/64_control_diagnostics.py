#!/usr/bin/env python3
"""Correct similarity and training-size diagnostics for the size-matched control.

The first version of the control report summarised FOLD MAXIMA by their mean and called the
result an achieved maximum. Three distinct quantities are computed and labelled separately here:

  A. max_over_test        -- the maximum, across ALL evaluated test complexes of all folds, of a
                             test complex's similarity to its fold's training set. This is the
                             only quantity entitled to be called "the achieved maximum".
  B. mean_of_fold_maxima  -- the mean over folds of each fold's maximum. Always <= A.
  C. mean_per_test_nn     -- the mean, over test complexes, of each complex's nearest-training
                             similarity. A different quantity again, and the smallest of the three.

Training sizes are reported per fold exactly; any average is labelled an average.

usage: python scripts/64_control_diagnostics.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from Levenshtein import ratio

sys.path.insert(0, ".")
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/"results/review/control_estimands"; OUT.mkdir(parents=True, exist_ok=True)
THRESHOLDS = [0.95, 0.85, 0.70]
DRAW_SEEDS = [1001, 1002, 1003, 1004, 1005]


def main():
    md = pd.read_csv(ROOT/"data/processed/complex_metadata.tsv", sep="\t")
    seqs = md.tcr_seq.tolist(); n = len(seqs)
    Sim = np.eye(n)
    for i in range(n):
        for j in range(i+1, n):
            Sim[i, j] = Sim[j, i] = ratio(seqs[i], seqs[j])
    folds = np.array_split(np.random.default_rng(0).permutation(n), 5)

    # reproduce the draws exactly as scripts/60 does
    k = {}
    for thr in THRESHOLDS:
        for fi, te in enumerate(folds):
            tr = np.setdiff1d(np.arange(n), te)
            k[(thr, fi)] = int(np.array([Sim[t, te].max() < thr for t in tr]).sum())
    subs = {}
    for ds in DRAW_SEEDS:
        rng = np.random.default_rng(ds)
        for thr in THRESHOLDS:
            for fi, te in enumerate(folds):
                tr = np.setdiff1d(np.arange(n), te)
                subs[(ds, thr, fi)] = np.sort(rng.choice(tr, k[(thr, fi)], replace=False))

    def summarise(train_of_fold):
        """train_of_fold: fold index -> training indices. Returns A, B, C."""
        fold_max, per_test = [], []
        for fi, te in enumerate(folds):
            tr = train_of_fold(fi)
            nn = np.array([Sim[i, tr].max() for i in te])   # each test complex's nearest training
            fold_max.append(nn.max()); per_test.append(nn)
        per_test = np.concatenate(per_test)
        return dict(max_over_test=float(per_test.max()),
                    mean_of_fold_maxima=float(np.mean(fold_max)),
                    mean_per_test_nn=float(per_test.mean()),
                    fold_maxima=[round(float(x), 4) for x in fold_max])

    rows = []
    for thr in THRESHOLDS:
        tr_sizes = [k[(thr, fi)] for fi in range(5)]
        filt = summarise(lambda fi: np.setdiff1d(np.arange(n), folds[fi])[
            np.array([Sim[t, folds[fi]].max() < thr
                      for t in np.setdiff1d(np.arange(n), folds[fi])])])
        rows.append(dict(regime="similarity-filtered", threshold=thr, draw_seed="-",
                         train_sizes_per_fold=",".join(map(str, tr_sizes)),
                         train_size_mean=round(float(np.mean(tr_sizes)), 1),
                         **{kk: (round(v, 4) if isinstance(v, float) else v)
                            for kk, v in filt.items()}))
        for ds in DRAW_SEEDS:
            s = summarise(lambda fi: subs[(ds, thr, fi)])
            rows.append(dict(regime="size-matched random", threshold=thr, draw_seed=ds,
                             train_sizes_per_fold=",".join(str(len(subs[(ds, thr, fi)])) for fi in range(5)),
                             train_size_mean=round(float(np.mean(tr_sizes)), 1),
                             **{kk: (round(v, 4) if isinstance(v, float) else v)
                                for kk, v in s.items()}))
    # unrestricted reference
    unr = summarise(lambda fi: np.setdiff1d(np.arange(n), folds[fi]))
    rows.append(dict(regime="unrestricted", threshold=1.01, draw_seed="-",
                     train_sizes_per_fold=",".join(str(n-len(f)) for f in folds),
                     train_size_mean=round(float(np.mean([n-len(f) for f in folds])), 1),
                     **{kk: (round(v, 4) if isinstance(v, float) else v) for kk, v in unr.items()}))

    t = pd.DataFrame(rows)
    t.to_csv(OUT/"similarity_diagnostics.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print("=== similarity diagnostics, three distinct quantities ===")
    print(t[["regime", "threshold", "draw_seed", "train_sizes_per_fold", "train_size_mean",
             "max_over_test", "mean_of_fold_maxima", "mean_per_test_nn"]].to_string(index=False))
    print("\n=== condensed, size-matched draws aggregated over the five draws ===")
    g = t[t.regime == "size-matched random"].groupby("threshold")
    print(pd.DataFrame(dict(
        max_over_test_overall=g.max_over_test.max(),
        mean_of_fold_maxima_across_draws=g.mean_of_fold_maxima.mean().round(4),
        mean_per_test_nn_across_draws=g.mean_per_test_nn.mean().round(4))).to_string())
    print(f"\nwrote {OUT/'similarity_diagnostics.tsv'}")


if __name__ == "__main__":
    main()
