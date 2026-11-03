#!/usr/bin/env python3
"""p82v1: TEIM-Res baseline, pairwise task, tcr70 only. Direct copy of
scripts/36_teim_res.py logic, pointed at the 82-position data.
"""
import json, sys, time
import numpy as np, pandas as pd, torch, torch.nn as nn
from pathlib import Path
sys.path.insert(0,'.')
from src.models.grip import encode_seq, SEG_ORDER, SEG_LENS, TCR_LEN
from src.eval.metrics import per_complex, cluster_bootstrap_diff
from src.eval.resume import refuse_overwrite

SEEDS=(42,43,44); EPOCHS=400
FINAL = "results/tables/teim_res_pairwise_p82v1.npz"

class SeqEncoder(nn.Module):
    def __init__(self, d=64, k=3, layers=2):
        super().__init__()
        from src.models.grip import VOCAB
        self.emb = nn.Embedding(VOCAB, d, padding_idx=0)
        blocks=[]
        for _ in range(layers):
            blocks += [nn.Conv1d(d, d, k, padding=k//2), nn.BatchNorm1d(d), nn.ReLU()]
        self.cnn = nn.Sequential(*blocks)
    def forward(self, x):
        h = self.emb(x).transpose(1,2)
        return self.cnn(h).transpose(1,2)

class ResBlock2D(nn.Module):
    def __init__(self, c, k=3, dil=1):
        super().__init__()
        self.c1=nn.Conv2d(c,c,k,padding=dil*(k//2),dilation=dil); self.b1=nn.BatchNorm2d(c)
        self.c2=nn.Conv2d(c,c,k,padding=dil*(k//2),dilation=dil); self.b2=nn.BatchNorm2d(c)
    def forward(self,x):
        h=torch.relu(self.b1(self.c1(x)))
        return torch.relu(x + self.b2(self.c2(h)))

class TEIMRes(nn.Module):
    def __init__(self, mhc_len, d=64, c=32):
        super().__init__()
        self.enc_t = SeqEncoder(d); self.enc_m = SeqEncoder(d)
        self.proj = nn.Conv2d(d, c, 1)
        self.blocks = nn.Sequential(ResBlock2D(c,3,1), ResBlock2D(c,3,2), ResBlock2D(c,3,4))
        self.out = nn.Conv2d(c, 1, 1)
    def forward(self, t, m, adj=None):
        ht = self.enc_t(t); hm = self.enc_m(m)
        inter = ht.unsqueeze(2) * hm.unsqueeze(1)
        x = self.proj(inter.permute(0,3,1,2))
        x = self.blocks(x)
        return self.out(x).squeeze(1)

def main():
    t0 = time.time()
    dev='cuda' if torch.cuda.is_available() else 'cpu'
    print(f"device: {dev}")
    refuse_overwrite(FINAL)
    Z=np.load("data/processed/pair_maps_p82v1.npz", allow_pickle=True)
    Y,S = Z["Y"].astype(np.float32), Z["S"].astype(np.float32)
    meta=pd.read_parquet("data/processed/pair_meta_p82v1.parquet")
    assert list(Z["pdb"]) == list(meta.pdb)
    md=pd.read_csv("data/processed/complex_metadata.tsv",sep="\t").set_index("pdb")
    meta=meta.join(md[["clu70"]],on="pdb")
    mhc_len=Y.shape[2]; assert mhc_len == 82
    T=np.stack([np.concatenate([encode_seq(r[k],SEG_LENS[k]) for k in SEG_ORDER]) for _,r in meta.iterrows()])
    M=np.stack([encode_seq(s,mhc_len) for s in meta.mhc_contact])
    proto=json.load(open("data/splits_v2/protocols.json"))["tcr70"]
    pdb2i={p:i for i,p in enumerate(meta.pdb)}

    P=np.full(Y.shape,np.nan,np.float32); used=np.zeros(len(Y),bool)
    for fold in proto:
        te=np.array([pdb2i[p] for p in fold if p in pdb2i])
        if len(te)<3: continue
        tr=np.setdiff1d(np.arange(len(Y)),te); used[te]=True
        assert set(tr.tolist()) & set(te.tolist()) == set()
        ps=[]
        for sd in SEEDS:
            torch.manual_seed(sd)
            net=TEIMRes(mhc_len).to(dev)
            Xa,Ma,Ya,Sa=(torch.tensor(z[tr]).to(dev) for z in (T,M,Y,S))
            pw=((Ya==0)&(Sa>0)).sum()/torch.clamp(((Ya==1)&(Sa>0)).sum(),min=1)
            opt=torch.optim.AdamW(net.parameters(),lr=1e-3,weight_decay=0.01)
            lf=nn.BCEWithLogitsLoss(reduction='none',pos_weight=pw)
            for ep in range(EPOCHS):
                net.train(); opt.zero_grad()
                l=(lf(net(Xa,Ma),Ya)*Sa).sum()/Sa.sum(); l.backward()
                nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
            net.eval()
            with torch.no_grad():
                ps.append(torch.sigmoid(net(torch.tensor(T[te]).to(dev),
                                            torch.tensor(M[te]).to(dev))).cpu().numpy())
        P[te]=np.mean(ps,0)
    idx=np.where(used)[0]
    mm=per_complex(Y[idx],P[idx],S[idx])
    print(f"  TEIM-Res  n={len(idx)}  macro-AUPRC={np.nanmean(mm['auprc']):.4f}  macro-AUC={np.nanmean(mm['auc']):.4f}")
    tmp = FINAL.replace(".npz", ".tmp.npz")   # np.savez auto-appends .npz if the path lacks it
    np.savez(tmp, pdb=meta.pdb.values[idx], clu70=meta.clu70.values[idx],
             auprc=mm["auprc"], auc=mm["auc"], pred=P[idx])
    Path(tmp).replace(FINAL)
    print(f"wrote {FINAL}  ({time.time()-t0:.1f}s)")

if __name__=="__main__":
    main()
