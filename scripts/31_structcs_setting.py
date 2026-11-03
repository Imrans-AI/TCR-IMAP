#!/usr/bin/env python3
"""GRIP vs Struct-TCS on the competitor's own STRUCT-CS complexes (n=88 shared)."""
import json, re, sys
import numpy as np, pandas as pd, torch, torch.nn as nn
from sklearn.metrics import roc_auc_score, average_precision_score
sys.path.insert(0,'.')
from src.models.grip import encode_seq, SEG_ORDER, SEG_LENS
from src.eval.metrics import per_complex, cluster_bootstrap_diff
_src = open("scripts/21_struct_tcs.py").read().split("def main()")[0].replace("import argparse, json, sys","")
exec(_src)

SEEDS=(42,43,44); EPOCHS=300
def fit(cls,kw,T,M,A,Y,S,tr,va,te,seed,dev):
    torch.manual_seed(seed); net=cls(**kw).to(dev)
    g=lambda i: tuple(torch.tensor(z[i]).to(dev) for z in (T,M,A,Y,S))
    Xa,Ma,Aa,Ya,Sa=g(tr); Xv,Mv,Av,Yv,Sv=g(va); Xe,Me,Ae,_,_=g(te)
    pw=((Ya==0)&(Sa>0)).sum()/torch.clamp(((Ya==1)&(Sa>0)).sum(),min=1)
    opt=torch.optim.Adam(net.parameters(),lr=1e-3,weight_decay=1e-4)
    lf=nn.BCEWithLogitsLoss(reduction='none',pos_weight=pw)
    best,bad,bs=-1,0,None
    for ep in range(EPOCHS):
        net.train(); opt.zero_grad()
        l=(lf(net(Xa,Ma,Aa),Ya)*Sa).sum()/Sa.sum(); l.backward()
        nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
        if ep%5==0:
            net.eval()
            with torch.no_grad(): pv=torch.sigmoid(net(Xv,Mv,Av))
            k=Sv.bool(); sc=average_precision_score(Yv[k].cpu(),pv[k].cpu())
            if sc>best: best,bad,bs=sc,0,{a:b.clone() for a,b in net.state_dict().items()}
            else:
                bad+=1
                if bad>10: break
    net.load_state_dict(bs); net.eval()
    with torch.no_grad(): return torch.sigmoid(net(Xe,Me,Ae)).cpu().numpy()

if __name__=="__main__":
    dev='cuda' if torch.cuda.is_available() else 'cpu'
    ids=set(re.findall(r'[0-9][a-zA-Z0-9]{3}', open("data/raw/sgtpmi/data/STRUCT-CS/PDBids.txt").read().lower()))
    import os
    _pm=os.environ.get("PAIRMAPS","data/processed/pair_maps.npz")
    Z=np.load(_pm); A=np.load("data/processed/tcr_graphs.npz")["A"]
    meta=pd.read_parquet("data/processed/pair_meta.parquet")
    md=pd.read_csv("data/processed/complex_metadata.tsv",sep="\t").set_index("pdb")
    meta=meta.join(md[["clu70"]],on="pdb")
    keep=meta.pdb.str.lower().isin(ids).values
    print(f"shared complexes: {keep.sum()}")
    Y=(Z["Y"].sum(-1)>0).astype(np.float32)[keep]; S=(Z["S"].sum(-1)>0).astype(np.float32)[keep]
    A=A[keep]; meta=meta[keep].reset_index(drop=True)
    mhc_len=Z["Y"].shape[2]
    T=np.stack([np.concatenate([encode_seq(r[k],SEG_LENS[k]) for k in SEG_ORDER]) for _,r in meta.iterrows()])
    M=np.stack([encode_seq(s,mhc_len) for s in meta.mhc_contact])
    MODELS={"Ours":(OursPerRes,dict(mhc_len=mhc_len)),
            "Struct-TCS":(TCS,dict(mode="struct",mhc_len=mhc_len)),
            "Seq-TCS":(TCS,dict(mode="seq",mhc_len=mhc_len))}
    rng=np.random.default_rng(42)
    for setting in ("their_random","our_tcr70"):
        n=len(Y)
        if setting=="their_random":
            perm=rng.permutation(n); folds=[perm[:int(n*0.25)]]     # single split, ~their ratio
        else:
            cl=meta.clu70.values; u=rng.permutation(pd.unique(cl))
            folds=[np.where(np.isin(cl,g))[0] for g in np.array_split(u,5)]
        P={k:np.full(Y.shape,np.nan,np.float32) for k in MODELS}
        used=np.zeros(n,bool)
        for te in folds:
            rest=np.setdiff1d(np.arange(n),te)
            vm=rng.random(len(rest))<0.15; tr,va=rest[~vm],rest[vm]
            if len(va)<3 or len(tr)<15: continue
            used[te]=True
            for name,(cls,kw) in MODELS.items():
                P[name][te]=np.mean([fit(cls,kw,T,M,A,Y,S,tr,va,te,s,dev) for s in SEEDS],0)
        i=np.where(used)[0]
        mc={k:per_complex(Y[i],P[k][i],S[i]) for k in MODELS}
        print(f"\n[{setting}] n_test={len(i)}")
        for k in MODELS: print(f"   {k:12s} macro-AUPRC={np.nanmean(mc[k]['auprc']):.4f}")
        for k in ("Struct-TCS","Seq-TCS"):
            d,lo,hi,p=cluster_bootstrap_diff(mc["Ours"]["auprc"],mc[k]["auprc"],meta.clu70.values[i])
            print(f"   Ours vs {k:12s} {d:+.4f} [{lo:+.4f},{hi:+.4f}] p={p:.4f}")
