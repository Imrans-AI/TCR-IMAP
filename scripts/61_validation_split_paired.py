#!/usr/bin/env python3
"""Part C: the validation-split comparison, reported as what it is.

Two things this produces, both ADDITIONAL analyses -- neither re-designates anything:

1. Achieved validation fractions and optimisation-training sizes under both procedures, per fold,
   for all protocols. The manuscript currently reports the cluster-disjoint fractions only; the
   comparison changes validation membership, validation fraction and optimisation-training size
   together, so all three are tabulated.

2. A DIRECT paired before/after difference. The existing tables give TCR-IMAP's macro-AUPRC under
   each procedure separately; they do not give the paired difference between the two procedures
   with its own uncertainty. The stored per-complex arrays allow it, so it is computed with the
   same paired cluster bootstrap used everywhere else (clu70, 20,000 resamples, seed 42).

Why this matters: a move from p = 0.001 to p = 0.088 in the primary comparison is NOT itself
evidence of a significant difference between the two procedures. That is a separate question and
needs the paired test computed here. The result is reported whichever way it comes out.

usage: python scripts/61_validation_split_paired.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, ".")
from src.eval.metrics import cluster_bootstrap_diff

ROOT = Path(__file__).resolve().parent.parent
VF = ROOT/"results/review/validation_fix"
OUT = ROOT/"results/review/v3_partC"; OUT.mkdir(parents=True, exist_ok=True)
N_BOOT, SEED = 20000, 42
PROTOCOLS = ["random", "double", "allele", "temporal", "tcr70"]
METHODS = ["Ours", "Struct-TCS", "Seq-TCS", "PosFreq", "NN"]
LABEL = {"Ours": "TCR-IMAP", "Struct-TCS": "Struct-TCS", "Seq-TCS": "Seq-TCS",
         "PosFreq": "position prior", "NN": "NN transfer"}
MIN_FOLD = 3


# ---- 1. achieved validation fractions and optimisation-training sizes ------------------------
def cluster_disjoint_validation(clu, rest, rng, frac=0.15):
    cl = clu[rest]; uc = np.unique(cl); target = frac*len(rest); sel, n = [], 0
    for c in uc[rng.permutation(len(uc))]:
        if n >= target: break
        sel.append(c); n += int((cl == c).sum())
    return rest[np.isin(cl, sel)]


def fractions():
    import json
    df = pd.read_csv(ROOT/"data/processed/complex_metadata.tsv", sep="\t")
    proto = json.load(open(ROOT/"data/splits_v2/protocols.json"))
    clu = df.clu70.values; pdb2i = {p: i for i, p in enumerate(df.pdb)}
    rows = []
    for pname in PROTOCOLS:
        rng_c = np.random.default_rng(0); rng_r = np.random.default_rng(0)
        for fi, fold in enumerate(proto[pname]):
            te = np.array([pdb2i[p] for p in fold if p in pdb2i])
            if len(te) < MIN_FOLD: continue
            rest = np.setdiff1d(np.arange(len(df)), te)
            # cluster-disjoint (reported configuration)
            va_c = cluster_disjoint_validation(clu, rest, rng_c, 0.15)
            tr_c = np.setdiff1d(rest, va_c)
            # random 15 % (superseded configuration), same generator discipline
            vm = rng_r.random(len(rest)) < 0.15
            tr_r, va_r = rest[~vm], rest[vm]
            rows.append(dict(
                protocol=pname, fold=fi, n_test=len(te), n_nontest=len(rest),
                cd_n_val=len(va_c), cd_val_pct=round(100*len(va_c)/len(rest), 1),
                cd_n_opt_train=len(tr_c), cd_val_clusters=len(set(clu[va_c])),
                cd_straddling_clusters=len(set(clu[tr_c]) & set(clu[va_c])),
                rnd_n_val=len(va_r), rnd_val_pct=round(100*len(va_r)/len(rest), 1),
                rnd_n_opt_train=len(tr_r), rnd_val_clusters=len(set(clu[va_r])),
                rnd_straddling_clusters=len(set(clu[tr_r]) & set(clu[va_r]))))
    t = pd.DataFrame(rows)
    t.to_csv(OUT/"validation_fractions.tsv", sep="\t", index=False)
    print("=== achieved validation fractions and optimisation-training sizes, per fold ===")
    print(t.to_string(index=False))
    print(f"\ncluster-disjoint: validation {t.cd_val_pct.min():.1f}-{t.cd_val_pct.max():.1f} % "
          f"(median {t.cd_val_pct.median():.1f} %), optimisation-training "
          f"{t.cd_n_opt_train.min()}-{t.cd_n_opt_train.max()} complexes, "
          f"straddling clusters {t.cd_straddling_clusters.min()}-{t.cd_straddling_clusters.max()}")
    print(f"random 15 %:      validation {t.rnd_val_pct.min():.1f}-{t.rnd_val_pct.max():.1f} % "
          f"(median {t.rnd_val_pct.median():.1f} %), optimisation-training "
          f"{t.rnd_n_opt_train.min()}-{t.rnd_n_opt_train.max()} complexes, "
          f"straddling clusters {t.rnd_straddling_clusters.min()}-{t.rnd_straddling_clusters.max()} "
          f"(median {t.rnd_straddling_clusters.median():.1f}) in {int((t.rnd_straddling_clusters>0).sum())} "
          f"of {len(t)} folds")
    return t


# ---- 2. direct paired before/after difference -------------------------------------------------
def paired():
    rows = []
    for pname in PROTOCOLS:
        fc = VF/f"clusterdisjoint_run1/percomplex_{pname}_p84v2.npz"
        fr = VF/f"randomsplit_deterministic/percomplex_{pname}_p84v2.npz"
        if not (fc.exists() and fr.exists()):
            print(f"  [{pname}] stored arrays missing, skipped"); continue
        C = np.load(fc, allow_pickle=True); R = np.load(fr, allow_pickle=True)
        assert list(C["pdb"]) == list(R["pdb"]), f"{pname}: complex order differs between runs"
        clu = C["clu70"]
        for m in METHODS:
            k = f"{m}_auprc"
            if k not in C.files or k not in R.files: continue
            a, b = C[k], R[k]          # cluster-disjoint minus random
            d, lo, hi, p = cluster_bootstrap_diff(a, b, clu, n_boot=N_BOOT, seed=SEED)
            rows.append(dict(protocol=pname, method=LABEL[m], n=len(a),
                             random_split=round(float(np.nanmean(b)), 4),
                             cluster_disjoint=round(float(np.nanmean(a)), 4),
                             paired_diff=round(float(d), 4),
                             lo=round(float(lo), 4), hi=round(float(hi), 4),
                             p=round(float(p), 4),
                             significant="yes" if p < 0.05 else "no"))
    t = pd.DataFrame(rows)
    t.to_csv(OUT/"paired_before_after.tsv", sep="\t", index=False)
    print("\n=== ADDITIONAL ANALYSIS: direct paired difference, cluster-disjoint minus random ===")
    print("(paired by complex, clu70 cluster bootstrap, 20,000 resamples, seed 42)")
    print(t.to_string(index=False))
    return t


if __name__ == "__main__":
    fractions(); paired()
    print(f"\nwrote {OUT}")
