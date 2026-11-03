#!/usr/bin/env python3
"""p82v1: ablations on tcr70 (primary protocol), pairwise task. Direct copy of
scripts/29_ablations_v2.py logic, pointed at the 82-position data, with
per-ablation resume-safe caching.
"""
import json, sys, time
import numpy as np, pandas as pd, torch, torch.nn as nn
from pathlib import Path
sys.path.insert(0,'.')
from src.models.grip import encode_seq, SEG_ORDER, SEG_LENS, TCR_LEN, GERMLINE_SEGS
from src.eval.metrics import per_complex, cluster_bootstrap_diff
from src.eval.resume import resumable_unit, refuse_overwrite

_src = open("scripts/19_train_pairmap.py").read().split("def main()")[0]
_src = _src.replace("import argparse, json, sys","").replace("class PairNet(","class PairNetCls(")
exec(_src)

SEEDS=(42,43,44); EPOCHS=400
CACHE_DIR = "results/tables/_partial/p82v1/29b_ablations"
FINAL_CSV = "results/tables/ablations_p82v1.csv"

def seg_slice(name):
    o=0
    for s in SEG_ORDER:
        if s==name: return slice(o,o+SEG_LENS[s])
        o+=SEG_LENS[s]

def apply_ablation(name, T, M, rng):
    T, M = T.copy(), M.copy()
    if name=="no_mhc":   M[:] = 1
    if name=="no_tcr":   T[T>0] = 1
    if name=="shuffle_tcr":
        for s in SEG_ORDER:
            sl = seg_slice(s)
            for i in range(len(T)):
                blk = T[i, sl]; nz = blk[blk>0]
                rng.shuffle(nz); blk[blk>0] = nz; T[i, sl] = blk
    if name=="shuffle_mhc":
        for i in range(len(M)): rng.shuffle(M[i])
    if name=="no_cdr3":
        for s in ("cdr3a","cdr3b"): T[:, seg_slice(s)] = 0
    if name=="no_germline":
        for s in GERMLINE_SEGS: T[:, seg_slice(s)] = 0
    return T, M

class NoCross(PairNetCls):
    def forward(self, t, m):
        x = self.tenc(self.emb(t) + self.seg(self.sid)[None] + self.pos)
        mm = self.menc(self.memb(m) + self.mpos)
        logits = torch.einsum('bid,bjd->bij', self.wt(x), self.wm(mm))/x.shape[-1]**0.5
        return logits*self.scale + self.pair_pos[None] + self.bias

class DotHead(PairNetCls):
    def forward(self, t, m):
        x = self.tenc(self.emb(t) + self.seg(self.sid)[None] + self.pos)
        mm = self.menc(self.memb(m) + self.mpos)
        a,_ = self.cross(x, mm, mm); x = self.norm(x+a)
        logits = torch.einsum('bid,bjd->bij', x, mm)/x.shape[-1]**0.5
        return logits*self.scale + self.pair_pos[None] + self.bias

def run(name, Y,S,T,M,meta,proto,dev):
    rng = np.random.default_rng(0)
    Ta, Ma = apply_ablation(name, T, M, rng)
    cls = NoCross if name=="no_cross" else DotHead if name=="dot_head" else PairNetCls
    pdb2i = {p:i for i,p in enumerate(meta.pdb)}
    P = np.full(Y.shape, np.nan, np.float32); used=np.zeros(len(Y),bool)
    for fold in proto:
        te = np.array([pdb2i[p] for p in fold if p in pdb2i])
        if len(te)<3: continue
        tr = np.setdiff1d(np.arange(len(Y)), te); used[te]=True
        assert set(tr.tolist()) & set(te.tolist()) == set()
        ps=[]
        for sd in SEEDS:
            torch.manual_seed(sd)
            net = cls(mhc_len=M.shape[1]).to(dev)
            Xa,Ma_,Ya,Sa = (torch.tensor(z[tr]).to(dev) for z in (Ta,Ma,Y,S))
            pw = ((Ya==0)&(Sa>0)).sum()/torch.clamp(((Ya==1)&(Sa>0)).sum(),min=1)
            opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=0.01)
            lf = nn.BCEWithLogitsLoss(reduction='none', pos_weight=pw)
            for ep in range(EPOCHS):
                net.train(); opt.zero_grad()
                l=(lf(net(Xa,Ma_),Ya)*Sa).sum()/Sa.sum(); l.backward()
                nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
            net.eval()
            with torch.no_grad():
                ps.append(torch.sigmoid(net(torch.tensor(Ta[te]).to(dev),
                                            torch.tensor(Ma[te]).to(dev))).cpu().numpy())
        P[te]=np.mean(ps,0)
    i=np.where(used)[0]
    return per_complex(Y[i],P[i],S[i]), i

def main():
    t_all = time.time()
    dev='cuda' if torch.cuda.is_available() else 'cpu'
    print(f"device: {dev}")
    refuse_overwrite(FINAL_CSV)
    Z=np.load("data/processed/pair_maps_p82v1.npz", allow_pickle=True)
    Y,S = Z["Y"].astype(np.float32), Z["S"].astype(np.float32)
    meta=pd.read_parquet("data/processed/pair_meta_p82v1.parquet")
    assert list(Z["pdb"]) == list(meta.pdb)
    md=pd.read_csv("data/processed/complex_metadata.tsv",sep="\t").set_index("pdb")
    meta=meta.join(md[["clu70"]],on="pdb")
    T=np.stack([np.concatenate([encode_seq(r[k],SEG_LENS[k]) for k in SEG_ORDER]) for _,r in meta.iterrows()])
    M=np.stack([encode_seq(s,Y.shape[2]) for s in meta.mhc_contact])
    proto=json.load(open("data/splits_v2/protocols.json"))["tcr70"]

    ABL = ("full","no_mhc","no_tcr","shuffle_tcr","shuffle_mhc","no_cdr3","no_germline","no_cross","dot_head")
    cached = {}
    for name in ABL:
        t0 = time.time()
        def compute():
            m, idx = run(name,Y,S,T,M,meta,proto,dev)
            return {"auprc": m["auprc"], "idx": idx}
        result, resumed = resumable_unit(CACHE_DIR, name, compute)
        cached[name] = result
        print(f"  {name:13s} macro-AUPRC={np.nanmean(result['auprc']):.4f} (resumed={resumed}, {time.time()-t0:.1f}s)")

    base = cached["full"]
    cl = meta.clu70.values[base["idx"]]
    rows=[dict(ablation="full", auprc=float(np.nanmean(base['auprc'])), diff=0.0, lo=0.0, hi=0.0, p=1.0)]
    for name in ABL[1:]:
        m = cached[name]
        d,lo,hi,p = cluster_bootstrap_diff(base["auprc"], m["auprc"], cl)
        print(f"  {name:13s} full-abl={d:+.4f} [{lo:+.4f},{hi:+.4f}] p={p:.4f}")
        rows.append(dict(ablation=name, auprc=float(np.nanmean(m['auprc'])), diff=d, lo=lo, hi=hi, p=p))
    df = pd.DataFrame(rows)
    tmp = FINAL_CSV + ".tmp"
    df.to_csv(tmp, index=False)
    Path(tmp).replace(FINAL_CSV)
    print(f"\nwrote {FINAL_CSV}  (total {time.time()-t_all:.1f}s)")

if __name__=="__main__":
    main()
