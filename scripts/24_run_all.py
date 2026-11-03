#!/usr/bin/env python3
"""Evaluate all models on the pairwise task across all five frozen protocols."""
import argparse, json, sys
import numpy as np, pandas as pd, torch, torch.nn as nn
sys.path.insert(0,'.')
from src.models.grip import encode_seq, SEG_ORDER, SEG_LENS, TCR_LEN
from src.eval.metrics import per_complex, cluster_bootstrap_diff
from importlib import import_module
_src = open("scripts/19_train_pairmap.py").read().split("def main()")[0]
_src = _src.replace("import argparse, json, sys", "").replace("class PairNet(", "class PairNetCls(")
exec(_src)

def load():
    Z = np.load("data/processed/pair_maps.npz")
    meta = pd.read_parquet("data/processed/pair_meta.parquet")
    md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
    meta = meta.join(md[["clu70","clu80","clu90","dup_family"]], on="pdb")
    T = np.stack([np.concatenate([encode_seq(r[k], SEG_LENS[k]) for k in SEG_ORDER]) for _,r in meta.iterrows()])
    M = np.stack([encode_seq(s, Z["Y"].shape[2]) for s in meta.mhc_contact])
    return Z["Y"].astype(np.float32), Z["S"].astype(np.float32), T, M, meta

def fit_predict(Y,S,T,M,tr,te,seed,epochs=400,dev='cuda'):
    torch.manual_seed(seed)
    net = PairNetCls(mhc_len=M.shape[1]).to(dev)
    Xa,Ma,Ya,Sa = (torch.tensor(z[tr]).to(dev) for z in (T,M,Y,S))
    pw = ((Ya==0)&(Sa>0)).sum()/torch.clamp(((Ya==1)&(Sa>0)).sum(),min=1)
    opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=0.01)
    lf = nn.BCEWithLogitsLoss(reduction='none', pos_weight=pw)
    for ep in range(epochs):
        net.train(); opt.zero_grad()
        l=(lf(net(Xa,Ma),Ya)*Sa).sum()/Sa.sum(); l.backward()
        nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
    net.eval()
    with torch.no_grad():
        return torch.sigmoid(net(torch.tensor(T[te]).to(dev), torch.tensor(M[te]).to(dev))).cpu().numpy()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--seeds",nargs="+",type=int,default=[42,43,44]); a=ap.parse_args()
    dev='cuda' if torch.cuda.is_available() else 'cpu'
    Y,S,T,M,meta = load()
    proto = json.load(open("data/splits_v2/protocols.json"))
    pdb2i = {p:i for i,p in enumerate(meta.pdb)}
    rows=[]
    for pname, folds in proto.items():
        preds = {k: np.full(Y.shape, np.nan, np.float32) for k in ("model","pair_freq")}
        covered = np.zeros(len(Y), bool)
        for fold in folds:
            te = np.array([pdb2i[p] for p in fold if p in pdb2i])
            if len(te)<3: continue
            tr = np.setdiff1d(np.arange(len(Y)), te)
            covered[te]=True
            f = (Y[tr]*S[tr]).sum(0)/np.maximum(S[tr].sum(0),1)
            preds["pair_freq"][te] = f[None]
            ps = [fit_predict(Y,S,T,M,tr,te,s,dev=dev) for s in a.seeds]
            preds["model"][te] = np.mean(ps,0)
        idx = np.where(covered)[0]
        clusters = meta.clu70.values[idx]
        np.savez_compressed(f"results/tables/pairmaps_{pname}.npz",
                            pdb=meta.pdb.values, clu70=meta.clu70.values,
                            covered=covered, Y=Y.astype(np.int8), S=S.astype(np.int8),
                            pred_model=preds["model"], pred_freq=preds["pair_freq"])
        print(f"  dumped results/tables/pairmaps_{pname}.npz")
        mm = {k: per_complex(Y[idx], preds[k][idx], S[idx]) for k in ("model","pair_freq")}
        d,lo,hi,p = cluster_bootstrap_diff(mm["model"]["auprc"], mm["pair_freq"]["auprc"], clusters)
        rows.append(dict(protocol=pname, n=len(idx),
            model_auprc=round(float(np.nanmean(mm["model"]["auprc"])),4),
            freq_auprc=round(float(np.nanmean(mm["pair_freq"]["auprc"])),4),
            model_p_at_L=round(float(np.nanmean(mm["model"]["p_at_L"])),4),
            freq_p_at_L=round(float(np.nanmean(mm["pair_freq"]["p_at_L"])),4),
            diff=round(d,4), ci_lo=round(lo,4), ci_hi=round(hi,4), pval=round(p,4)))
        print(f"{pname:9s} n={len(idx):3d} macro-AUPRC model={rows[-1]['model_auprc']:.4f} "
              f"freq={rows[-1]['freq_auprc']:.4f} diff={d:+.4f} [{lo:+.4f},{hi:+.4f}] p={p:.4f}")
    pd.DataFrame(rows).to_csv("results/tables/protocol_main.csv", index=False)

if __name__=="__main__":
    main()
