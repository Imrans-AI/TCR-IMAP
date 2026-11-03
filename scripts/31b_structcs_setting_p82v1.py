#!/usr/bin/env python3
"""p82v1: GRIP vs Struct-TCS/Seq-TCS on the comparator's dataset subset (88 of their 94
STRUCT-CS complexes shared with our 174-complex corpus), at the comparator's reported
cutoff. See results/review/PRIMARY_82_POSITION_RERUN_PLAN.md section 13 for exactly which
dimensions are and are not actually matched to the original comparator evaluation -- only
the dataset subset and the contact cutoff are; the MHC candidate-position set, TCR
representation, task pooling, and split scheme remain ours throughout, in both the
28-position and this 82-position analysis.

Unlike the original scripts/31_structcs_setting.py (which only printed to stdout), this
version writes a results table, closing the gap noted in the main audit (section 9.3).
"""
import json, re, sys, time, os
import numpy as np, pandas as pd, torch, torch.nn as nn
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score
sys.path.insert(0,'.')
from src.models.grip import encode_seq, SEG_ORDER, SEG_LENS
from src.eval.metrics import per_complex, cluster_bootstrap_diff
from src.eval.resume import resumable_unit, refuse_overwrite

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

def load(pairmaps_path):
    ids=set(re.findall(r'[0-9][a-zA-Z0-9]{3}', open("data/raw/sgtpmi/data/STRUCT-CS/PDBids.txt").read().lower()))
    Z=np.load(pairmaps_path, allow_pickle=True); A=np.load("data/processed/tcr_graphs.npz")["A"]
    meta_path = pairmaps_path.replace("pair_maps", "pair_meta").replace(".npz", ".parquet")
    meta=pd.read_parquet(meta_path)
    assert list(Z["pdb"]) == list(meta.pdb)
    md=pd.read_csv("data/processed/complex_metadata.tsv",sep="\t").set_index("pdb")
    meta=meta.join(md[["clu70"]],on="pdb")
    keep=meta.pdb.str.lower().isin(ids).values
    print(f"shared complexes: {keep.sum()}")
    Y=(Z["Y"].sum(-1)>0).astype(np.float32)[keep]; S=(Z["S"].sum(-1)>0).astype(np.float32)[keep]
    A=A[keep]; meta=meta[keep].reset_index(drop=True)
    mhc_len=Z["Y"].shape[2]
    T=np.stack([np.concatenate([encode_seq(r[k],SEG_LENS[k]) for k in SEG_ORDER]) for _,r in meta.iterrows()])
    M=np.stack([encode_seq(s,mhc_len) for s in meta.mhc_contact])
    return T,M,A,Y,S,meta,keep.sum()

def run_setting(setting, T,M,A,Y,S,meta,dev):
    MODELS={"Ours":(OursPerRes,dict(mhc_len=M.shape[1])),
            "Struct-TCS":(TCS,dict(mode="struct",mhc_len=M.shape[1])),
            "Seq-TCS":(TCS,dict(mode="seq",mhc_len=M.shape[1]))}
    rng=np.random.default_rng(42)
    n=len(Y)
    if setting=="their_random":
        perm=rng.permutation(n); folds=[perm[:int(n*0.25)]]
    else:
        cl=meta.clu70.values; u=rng.permutation(pd.unique(cl))
        folds=[np.where(np.isin(cl,g))[0] for g in np.array_split(u,5)]
    P={k:np.full(Y.shape,np.nan,np.float32) for k in MODELS}
    used=np.zeros(n,bool)
    for te in folds:
        rest=np.setdiff1d(np.arange(n),te)
        vm=rng.random(len(rest))<0.15; tr,va=rest[~vm],rest[vm]
        if len(va)<3 or len(tr)<15: continue
        assert not (np.isin(np.arange(n),te) & np.isin(np.arange(n),tr)).any()
        used[te]=True
        for name,(cls,kw) in MODELS.items():
            P[name][te]=np.mean([fit(cls,kw,T,M,A,Y,S,tr,va,te,s,dev) for s in SEEDS],0)
    i=np.where(used)[0]
    mc={k:per_complex(Y[i],P[k][i],S[i]) for k in MODELS}
    line = dict(setting=setting, n_test=len(i))
    for k in MODELS: line[f"{k}_auprc"] = round(float(np.nanmean(mc[k]["auprc"])),4)
    for k in ("Struct-TCS","Seq-TCS"):
        d,lo,hi,p=cluster_bootstrap_diff(mc["Ours"]["auprc"],mc[k]["auprc"],meta.clu70.values[i])
        line[f"vs_{k}"] = f"{d:+.4f} [{lo:+.4f},{hi:+.4f}] p={p:.4f}"
    return line

def run_cutoff(tag, pairmaps_path, cache_dir, dev):
    T,M,A,Y,S,meta,n_shared = load(pairmaps_path)
    rows=[]
    for setting in ("their_random","our_tcr70"):
        t0=time.time()
        def compute():
            line = run_setting(setting, T,M,A,Y,S,meta,dev)
            return {"__json__": np.array(json.dumps(line))}
        result, resumed = resumable_unit(cache_dir, f"{tag}_{setting}", compute)
        line = json.loads(str(result["__json__"]))
        line["cutoff_tag"] = tag; line["n_shared_with_comparator"] = int(n_shared)
        rows.append(line)
        print(f"[{tag} / {setting}] n_test={line['n_test']} Ours={line['Ours_auprc']:.4f} "
              f"(resumed={resumed}, {time.time()-t0:.1f}s)")
    return rows

def main():
    t_all=time.time()
    dev='cuda' if torch.cuda.is_available() else 'cpu'
    print(f"device: {dev}")
    out_45 = "results/tables/structcs_setting_p82v1_4p5A.csv"
    out_6A = "results/tables/structcs_setting_p82v1_6A.csv"
    refuse_overwrite(out_45); refuse_overwrite(out_6A)

    rows_45 = run_cutoff("4p5A", "data/processed/pair_maps_p82v1.npz",
                          "results/tables/_partial/p82v1/31b_structcs_45A", dev)
    df45 = pd.DataFrame(rows_45)
    tmp = out_45 + ".tmp"; df45.to_csv(tmp, index=False); Path(tmp).replace(out_45)
    print(f"wrote {out_45}")

    rows_6A = run_cutoff("6A", "data/processed/pair_maps_p82v1_6A.npz",
                         "results/tables/_partial/p82v1/31b_structcs_6A", dev)
    df6A = pd.DataFrame(rows_6A)
    tmp = out_6A + ".tmp"; df6A.to_csv(tmp, index=False); Path(tmp).replace(out_6A)
    print(f"wrote {out_6A}")
    print(f"\ntotal time: {time.time()-t_all:.1f}s")

if __name__=="__main__":
    main()
