#!/usr/bin/env python3
"""p82v1: per-residue task benchmark under all five frozen protocols.
Ours / Struct-TCS / Seq-TCS / position-frequency / nearest-neighbour transfer.
Direct copy of scripts/28_benchmark_v2.py logic, pointed at the 82-position data,
with per-protocol resume-safe caching (write-temp -> validate -> atomic rename).
"""
import json, sys, time
import numpy as np, pandas as pd, torch, torch.nn as nn
from pathlib import Path
from Levenshtein import ratio
sys.path.insert(0,'.')
from src.models.grip import encode_seq, SEG_ORDER, SEG_LENS, TCR_LEN
from src.eval.metrics import per_complex, cluster_bootstrap_diff
from src.eval.resume import resumable_unit, refuse_overwrite

_src = open("scripts/21_struct_tcs.py").read().split("def main()")[0]
_src = _src.replace("import argparse, json, sys","")
exec(_src)   # gives GatedGAT, AttnBlock, GATBlock, TCS, OursPerRes

SEEDS = (42,43,44)
EPOCHS = 300
CACHE_DIR = "results/tables/_partial/p82v1/28b_benchmark"
FINAL_CSV = "results/tables/benchmark_p82v1.csv"

def fit(model_cls, kw, T,M,A,Y,S, tr, va, te, seed, dev):
    torch.manual_seed(seed)
    net = model_cls(**kw).to(dev)
    g = lambda i: tuple(torch.tensor(z[i]).to(dev) for z in (T,M,A,Y,S))
    Xa,Ma,Aa,Ya,Sa = g(tr); Xv,Mv,Av,Yv,Sv = g(va); Xe,Me,Ae,_,_ = g(te)
    pw = ((Ya==0)&(Sa>0)).sum()/torch.clamp(((Ya==1)&(Sa>0)).sum(),min=1)
    opt = torch.optim.Adam(net.parameters(), lr=1e-3, weight_decay=1e-4)
    lf = nn.BCEWithLogitsLoss(reduction='none', pos_weight=pw)
    from sklearn.metrics import average_precision_score
    best,bad,bs = -1,0,None
    for ep in range(EPOCHS):
        net.train(); opt.zero_grad()
        l=(lf(net(Xa,Ma,Aa),Ya)*Sa).sum()/Sa.sum(); l.backward()
        nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
        if ep%5==0:
            net.eval()
            with torch.no_grad(): pv=torch.sigmoid(net(Xv,Mv,Av))
            k=Sv.bool()
            sc=average_precision_score(Yv[k].cpu(),pv[k].cpu())
            if sc>best: best,bad,bs=sc,0,{a:b.clone() for a,b in net.state_dict().items()}
            else:
                bad+=1
                if bad>10: break
    net.load_state_dict(bs); net.eval()
    with torch.no_grad(): return torch.sigmoid(net(Xe,Me,Ae)).cpu().numpy()

def load():
    Z = np.load("data/processed/pair_maps_p82v1.npz", allow_pickle=True)
    A = np.load("data/processed/tcr_graphs.npz")["A"]
    meta = pd.read_parquet("data/processed/pair_meta_p82v1.parquet")
    assert list(Z["pdb"]) == list(meta.pdb), "pdb order mismatch, npz vs parquet"
    md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
    meta = meta.join(md[["clu70","tcr_seq"]], on="pdb")
    Y = (Z["Y"].sum(-1)>0).astype(np.float32)
    S = (Z["S"].sum(-1)>0).astype(np.float32)
    mhc_len = Z["Y"].shape[2]
    assert mhc_len == 82
    T = np.stack([np.concatenate([encode_seq(r[k], SEG_LENS[k]) for k in SEG_ORDER]) for _,r in meta.iterrows()])
    M = np.stack([encode_seq(s, mhc_len) for s in meta.mhc_contact])
    return Y, S, T, M, A, meta, mhc_len

def run_protocol(pname, folds, Y,S,T,M,A,meta,mhc_len,dev):
    MODELS = {"Ours": (OursPerRes, dict(mhc_len=mhc_len)),
              "Struct-TCS": (TCS, dict(mode="struct", mhc_len=mhc_len)),
              "Seq-TCS": (TCS, dict(mode="seq", mhc_len=mhc_len))}
    P = {k: np.full(Y.shape, np.nan, np.float32) for k in list(MODELS)+["PosFreq","NN"]}
    used = np.zeros(len(Y), bool)
    rng = np.random.default_rng(0)
    for fold in folds:
        pdb2i = {p:i for i,p in enumerate(meta.pdb)}
        te = np.array([pdb2i[p] for p in fold if p in pdb2i])
        if len(te)<3: continue
        rest = np.setdiff1d(np.arange(len(Y)), te)
        vmask = rng.random(len(rest))<0.15
        tr, va = rest[~vmask], rest[vmask]
        if len(va)<3 or len(tr)<20: continue
        assert set(tr.tolist()) & set(te.tolist()) == set()
        assert set(va.tolist()) & set(te.tolist()) == set()
        used[te]=True
        P["PosFreq"][te] = ((Y[rest]*S[rest]).sum(0)/np.maximum(S[rest].sum(0),1))[None]
        for i in te:
            j = rest[int(np.argmax([ratio(meta.tcr_seq.iloc[i], meta.tcr_seq.iloc[k]) for k in rest]))]
            P["NN"][i] = Y[j]
        for name,(cls,kw) in MODELS.items():
            ps=[fit(cls,kw,T,M,A,Y,S,tr,va,te,s,dev) for s in SEEDS]
            P[name][te]=np.mean(ps,0)
    idx=np.where(used)[0]; cl=meta.clu70.values[idx]
    mc={k: per_complex(Y[idx],P[k][idx],S[idx]) for k in P}
    line={"protocol":pname,"n":len(idx)}
    for k in P: line[k]=round(float(np.nanmean(mc[k]["auprc"])),4)
    for k in P:
        if k!="Ours":
            d,lo,hi,p = cluster_bootstrap_diff(mc["Ours"]["auprc"], mc[k]["auprc"], cl)
            line[f"vs_{k}"]=f"{d:+.4f} [{lo:+.4f},{hi:+.4f}] p={p:.4f}"
    percomplex = dict(pdb=meta.pdb.values[idx], clu70=meta.clu70.values[idx],
                      allele=meta.allele.values[idx].astype(str),
                      **{f"{k}_auprc": mc[k]["auprc"] for k in P},
                      **{f"{k}_auc": mc[k]["auc"] for k in P})
    return line, percomplex

def main():
    dev='cuda' if torch.cuda.is_available() else 'cpu'
    print(f"device: {dev}")
    refuse_overwrite(FINAL_CSV)
    Y, S, T, M, A, meta, mhc_len = load()
    proto = json.load(open("data/splits_v2/protocols.json"))
    rows = []
    for pname, folds in proto.items():
        t0 = time.time()
        def compute():
            line, percomplex = run_protocol(pname, folds, Y,S,T,M,A,meta,mhc_len,dev)
            out = dict(percomplex)
            for k,v in line.items():
                out[f"__line__{k}"] = v
            return out
        def validate(d):
            # per-complex arrays here are already restricted to the covered/used subset
            # (see run_protocol's `idx = np.where(used)[0]` before building `percomplex`),
            # so a plain all-finite check is correct here -- unlike 24b's full-length
            # pred_model/pred_freq arrays, which legitimately contain NaN for uncovered rows.
            for k, v in d.items():
                if k.startswith("__line__"): continue
                a = np.asarray(v)
                if a.dtype.kind == "f" and not np.isfinite(a).all():
                    return False
            return True
        result, resumed = resumable_unit(CACHE_DIR, f"protocol_{pname}", compute, validate_fn=validate)
        dt = time.time()-t0
        print(f"[{pname}] n={int(result['__line__n'])} Ours_auprc={result['__line__Ours']:.4f} "
              f"(resumed={resumed}, {dt:.1f}s)")
        line = {k[len('__line__'):]: v for k,v in result.items() if k.startswith('__line__')}
        for k in ('n',): line[k] = int(line[k])
        rows.append(line)
        percomplex = {k: v for k,v in result.items() if not k.startswith('__line__')}
        pc_path = f"results/tables/percomplex_{pname}_p82v1.npz"
        if Path(pc_path).exists():
            print(f"  {pc_path} already finalized, skipping rewrite (resume)")
        else:
            pc_tmp = pc_path.replace(".npz", ".tmp.npz")
            np.savez_compressed(pc_tmp, **percomplex)
            Path(pc_tmp).replace(pc_path)
            print(f"  wrote {pc_path}")

    df = pd.DataFrame(rows)
    tmp = FINAL_CSV + ".tmp"
    df.to_csv(tmp, index=False)
    Path(tmp).replace(FINAL_CSV)
    print(f"\nwrote {FINAL_CSV}")
    print(df.to_string(index=False))

if __name__=="__main__":
    main()
