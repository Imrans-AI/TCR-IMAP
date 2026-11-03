#!/usr/bin/env python3
"""
Distance-to-training analysis.

For a sweep of maximum allowed train/test TCR identity thresholds, build splits where
no test complex has a training neighbour above the threshold, then measure model and
frequency-prior performance. Expectation: the model's advantage decays with the
threshold while the frequency prior stays flat.
"""
import json, sys
import numpy as np, pandas as pd, torch, torch.nn as nn
from Levenshtein import ratio
sys.path.insert(0,'.')
from src.models.grip import encode_seq, SEG_ORDER, SEG_LENS
from src.eval.metrics import per_complex, cluster_bootstrap_diff

_src = open("scripts/19_train_pairmap.py").read().split("def main()")[0]
_src = _src.replace("import argparse, json, sys","").replace("class PairNet(","class PairNetCls(")
exec(_src)

THRESHOLDS = [1.01, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70]
SEEDS = [42, 43, 44]

def main():
    dev='cuda' if torch.cuda.is_available() else 'cpu'
    Z = np.load("data/processed/pair_maps.npz")
    Y, S = Z["Y"].astype(np.float32), Z["S"].astype(np.float32)
    meta = pd.read_parquet("data/processed/pair_meta.parquet")
    md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
    meta = meta.join(md[["clu70","tcr_seq"]], on="pdb")
    T = np.stack([np.concatenate([encode_seq(r[k], SEG_LENS[k]) for k in SEG_ORDER]) for _,r in meta.iterrows()])
    M = np.stack([encode_seq(s, Y.shape[2]) for s in meta.mhc_contact])
    seqs = meta.tcr_seq.tolist(); n = len(seqs)

    print("computing pairwise identity matrix ...", flush=True)
    Sim = np.eye(n)
    for i in range(n):
        for j in range(i+1, n):
            Sim[i,j] = Sim[j,i] = ratio(seqs[i], seqs[j])

    rng = np.random.default_rng(0)
    folds = np.array_split(rng.permutation(n), 5)
    rows = []
    for thr in THRESHOLDS:
        preds_m = np.full(Y.shape, np.nan, np.float32)
        preds_f = np.full(Y.shape, np.nan, np.float32)
        used = np.zeros(n, bool)
        for te in folds:
            tr = np.setdiff1d(np.arange(n), te)
            # drop training complexes too similar to any test complex
            keep = np.array([Sim[t, te].max() < thr for t in tr])
            tr2 = tr[keep]
            if len(tr2) < 30: continue
            used[te] = True
            f = (Y[tr2]*S[tr2]).sum(0)/np.maximum(S[tr2].sum(0),1)
            preds_f[te] = f[None]
            ps = []
            for sd in SEEDS:
                torch.manual_seed(sd)
                net = PairNetCls(mhc_len=M.shape[1]).to(dev)
                Xa,Ma,Ya,Sa = (torch.tensor(z[tr2]).to(dev) for z in (T,M,Y,S))
                pw = ((Ya==0)&(Sa>0)).sum()/torch.clamp(((Ya==1)&(Sa>0)).sum(),min=1)
                opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=0.01)
                lf = nn.BCEWithLogitsLoss(reduction='none', pos_weight=pw)
                for ep in range(400):
                    net.train(); opt.zero_grad()
                    l=(lf(net(Xa,Ma),Ya)*Sa).sum()/Sa.sum(); l.backward()
                    nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
                net.eval()
                with torch.no_grad():
                    ps.append(torch.sigmoid(net(torch.tensor(T[te]).to(dev),
                                                torch.tensor(M[te]).to(dev))).cpu().numpy())
            preds_m[te] = np.mean(ps,0)
        idx = np.where(used)[0]
        mm = per_complex(Y[idx], preds_m[idx], S[idx])
        ff = per_complex(Y[idx], preds_f[idx], S[idx])
        d,lo,hi,p = cluster_bootstrap_diff(mm["auprc"], ff["auprc"], meta.clu70.values[idx])
        # mean realised max train-test similarity at this threshold
        realised = np.mean([Sim[i, np.setdiff1d(np.arange(n),[i])][Sim[i, np.setdiff1d(np.arange(n),[i])] < thr].max()
                            if (Sim[i, np.setdiff1d(np.arange(n),[i])] < thr).any() else 0 for i in idx])
        rows.append(dict(threshold=thr, n=len(idx),
                         model=round(float(np.nanmean(mm["auprc"])),4),
                         freq=round(float(np.nanmean(ff["auprc"])),4),
                         diff=round(d,4), lo=round(lo,4), hi=round(hi,4), p=round(p,4),
                         mean_max_sim=round(float(realised),3)))
        print(f"  thr<{thr:<5} n={len(idx):3d} model={rows[-1]['model']:.4f} freq={rows[-1]['freq']:.4f} "
              f"diff={d:+.4f} [{lo:+.4f},{hi:+.4f}] p={p:.4f}", flush=True)
    pd.DataFrame(rows).to_csv("results/tables/similarity_curve.csv", index=False)
    print("\nwrote results/tables/similarity_curve.csv")

if __name__=="__main__":
    main()
