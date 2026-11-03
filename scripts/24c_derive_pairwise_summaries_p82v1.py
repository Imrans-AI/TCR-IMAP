#!/usr/bin/env python3
"""p82v1: derive the pairwise case-study summary tables (pairwise_ours_p82v1.npz,
percomplex_tcr70_paired_p82v1.tsv) from results/tables/pairmaps_tcr70_p82v1.npz
(produced by scripts/24b_run_all_p82v1.py). No training performed here -- this is a
pure post-processing/aggregation step.

Closes the reproducibility gap noted in the main audit (section 9.3): the 28-position
equivalents of these two files (pairwise_ours.npz, percomplex_tcr70_paired.tsv) are read
by several existing scripts but written by none of them -- no committed script produces
them. This script is the explicit, committed derivation for the p82v1 versions.
"""
import sys, time
import numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0,'.')
from src.eval.metrics import per_complex
from src.eval.resume import refuse_overwrite

FINAL_OURS = "results/tables/pairwise_ours_p82v1.npz"
FINAL_PAIRED = "results/tables/percomplex_tcr70_paired_p82v1.tsv"

def main():
    t0 = time.time()
    refuse_overwrite(FINAL_OURS); refuse_overwrite(FINAL_PAIRED)
    pm = np.load("results/tables/pairmaps_tcr70_p82v1.npz", allow_pickle=True)
    Y, S, pred_model, pred_freq = pm["Y"].astype(np.float32), pm["S"].astype(np.float32), pm["pred_model"], pm["pred_freq"]
    covered = pm["covered"]
    idx = np.where(covered)[0]
    pdbs = pm["pdb"][idx]; clu70 = pm["clu70"][idx]

    mm_model = per_complex(Y[idx], pred_model[idx], S[idx])
    mm_freq = per_complex(Y[idx], pred_freq[idx], S[idx])

    tmp1 = FINAL_OURS.replace(".npz", ".tmp.npz")
    np.savez(tmp1, pdb=pdbs, clu70=clu70, auprc=mm_model["auprc"], auc=mm_model["auc"])
    Path(tmp1).replace(FINAL_OURS)
    print(f"wrote {FINAL_OURS}  (n={len(idx)})")

    meta = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
    rows = []
    for k, pdb in enumerate(pdbs):
        m = meta.loc[pdb]
        rows.append(dict(pdb=pdb, clu70=int(clu70[k]),
                          grip=float(mm_model["auprc"][k]), prior=float(mm_freq["auprc"][k]),
                          delta=float(mm_model["auprc"][k]-mm_freq["auprc"][k]),
                          npos=int(mm_model["n_pos"][k]),
                          rate=float(mm_model["n_pos"][k])/(S[idx][k].sum()) if S[idx][k].sum()>0 else float("nan"),
                          allele=m.allele, peptide=m.peptide, trav=m.trav, trbv=m.trbv))
    df = pd.DataFrame(rows)
    assert df.pdb.nunique() == len(df) == len(idx)
    tmp2 = FINAL_PAIRED + ".tmp"
    df.to_csv(tmp2, sep="\t", index=False)
    Path(tmp2).replace(FINAL_PAIRED)
    print(f"wrote {FINAL_PAIRED}  (n={len(df)})")

    win = (df.delta>0).mean()
    print(f"\nGRIP > prior on {int((df.delta>0).sum())}/{len(df)} complexes ({100*win:.1f}%)")
    top = df.sort_values("delta", ascending=False)
    print("\ntop 3 wins (largest positive delta):")
    print(top.head(3)[["pdb","grip","prior","delta"]].to_string(index=False))
    print("\ntop 3 losses (largest negative delta):")
    print(top.tail(3)[["pdb","grip","prior","delta"]].to_string(index=False))
    print(f"\ndone in {time.time()-t0:.1f}s")

if __name__=="__main__":
    main()
