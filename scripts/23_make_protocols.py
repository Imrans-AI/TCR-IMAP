#!/usr/bin/env python3
"""Freeze five evaluation protocols + train/test similarity audit. Run once."""
import sys, json
import numpy as np, pandas as pd
from Levenshtein import ratio
from pathlib import Path

SEED = 42
CUTOFF = "2020-01-01"

def kfold_groups(groups, k, rng):
    u = pd.unique(groups); rng.shuffle(u)
    return [set(x) for x in np.array_split(u, k)]

def main():
    df = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t")
    df["dep_date"] = pd.to_datetime(df.dep_date)
    rng = np.random.default_rng(SEED)
    out = Path("data/splits_v2"); out.mkdir(exist_ok=True)
    proto = {}

    # 1 random (complex-level)
    r = rng.permutation(len(df))
    proto["random"] = [set(df.pdb.values[x]) for x in np.array_split(r, 5)]
    # 2 allele holdout
    proto["allele"] = [set(df.pdb[df.allele.isin(g)]) for g in kfold_groups(df.allele.fillna("NA").values, 5, rng)]
    # 3 TCR-cluster holdout at 70%
    proto["tcr70"] = [set(df.pdb[df.clu70.isin(g)]) for g in kfold_groups(df.clu70.values, 5, rng)]
    # 4 double holdout: cluster AND allele both unseen
    # double holdout: 3 folds so each is usable; test = unseen cluster AND unseen allele
    folds = []
    gcs = kfold_groups(df.clu70.values,3,rng); gas = kfold_groups(df.allele.fillna("NA").values,3,rng)
    for gc, ga in zip(gcs, gas):
        te = df[df.clu70.isin(gc) & df.allele.isin(ga)]
        folds.append(set(te.pdb))
    proto["double"] = folds
    # 5 temporal (single split, not CV)
    proto["temporal"] = [set(df.pdb[df.dep_date >= CUTOFF])]

    audit = []
    for name, folds in proto.items():
        for fi, te in enumerate(folds):
            tr = df[~df.pdb.isin(te)]; ted = df[df.pdb.isin(te)]
            if len(ted)==0: continue
            # max TCR identity from each test complex to any training complex
            sims = [max((ratio(a,b) for b in tr.tcr_seq), default=0) for a in ted.tcr_seq]
            audit.append(dict(protocol=name, fold=fi, n_train=len(tr), n_test=len(ted),
                              max_sim_mean=round(float(np.mean(sims)),3),
                              max_sim_median=round(float(np.median(sims)),3),
                              frac_over_90=round(float(np.mean(np.array(sims)>0.9)),3),
                              test_alleles_unseen=int(ted.allele.nunique()-len(set(ted.allele)&set(tr.allele))),
                              dup_family_leak=int(len(set(ted.dup_family)&set(tr.dup_family)))))
    a = pd.DataFrame(audit)
    print(a.groupby("protocol")[["n_test","max_sim_mean","frac_over_90","dup_family_leak"]].mean().round(3).to_string())
    a.to_csv("data/splits_v2/similarity_audit.tsv", sep="\t", index=False)
    json.dump({k:[sorted(f) for f in v] for k,v in proto.items()},
              open("data/splits_v2/protocols.json","w"), indent=1)
    print("\nwrote data/splits_v2/protocols.json + similarity_audit.tsv")

if __name__=="__main__":
    main()
