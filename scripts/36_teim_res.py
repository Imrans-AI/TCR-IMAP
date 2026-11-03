#!/usr/bin/env python3
"""
TEIM-Res baseline, reimplemented for CDR-MHC pairwise contact maps.

Original (Peng et al., Nat Mach Intell 2023) predicts residue-level CDR3b-peptide
contact maps: per-sequence 1D CNN encoders -> outer-product interaction tensor ->
2D CNN with residual blocks -> per-pair contact score. We keep that architecture and
substitute our inputs (84 CDR positions x N MHC contact positions). Pretrained
peptide-autoencoder weights are not applicable here, so encoders are trained from
scratch; this is stated in Methods.

Evaluated on tcr70 with per-complex macro-AUPRC and cluster bootstrap, identically
to our model.
"""
import json, sys
import numpy as np, pandas as pd, torch, torch.nn as nn
sys.path.insert(0,'.')
from src.models.grip import encode_seq, VOCAB, SEG_ORDER, SEG_LENS, TCR_LEN
from src.eval.metrics import per_complex, cluster_bootstrap_diff

SEEDS=(42,43,44); EPOCHS=400

class SeqEncoder(nn.Module):
    """1D CNN sequence encoder (TEIM style)."""
    def __init__(self, d=64, k=3, layers=2):
        super().__init__()
        self.emb = nn.Embedding(VOCAB, d, padding_idx=0)
        blocks=[]
        for _ in range(layers):
            blocks += [nn.Conv1d(d, d, k, padding=k//2), nn.BatchNorm1d(d), nn.ReLU()]
        self.cnn = nn.Sequential(*blocks)
    def forward(self, x):
        h = self.emb(x).transpose(1,2)
        return self.cnn(h).transpose(1,2)            # B, L, d

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
        self.proj = nn.Conv2d(d, c, 1)               # from outer-product channels
        self.blocks = nn.Sequential(ResBlock2D(c,3,1), ResBlock2D(c,3,2), ResBlock2D(c,3,4))
        self.out = nn.Conv2d(c, 1, 1)
    def forward(self, t, m, adj=None):
        ht = self.enc_t(t)                            # B, T, d
        hm = self.enc_m(m)                            # B, M, d
        inter = ht.unsqueeze(2) * hm.unsqueeze(1)     # B, T, M, d  (outer product)
        x = self.proj(inter.permute(0,3,1,2))         # B, c, T, M
        x = self.blocks(x)
        return self.out(x).squeeze(1)                 # B, T, M

def main():
    dev='cuda' if torch.cuda.is_available() else 'cpu'
    Z=np.load("data/processed/pair_maps.npz")
    Y,S = Z["Y"].astype(np.float32), Z["S"].astype(np.float32)
    meta=pd.read_parquet("data/processed/pair_meta.parquet")
    md=pd.read_csv("data/processed/complex_metadata.tsv",sep="\t").set_index("pdb")
    meta=meta.join(md[["clu70"]],on="pdb")
    mhc_len=Y.shape[2]
    T=np.stack([np.concatenate([encode_seq(r[k],SEG_LENS[k]) for k in SEG_ORDER]) for _,r in meta.iterrows()])
    M=np.stack([encode_seq(s,mhc_len) for s in meta.mhc_contact])
    proto=json.load(open("data/splits_v2/protocols.json"))["tcr70"]
    pdb2i={p:i for i,p in enumerate(meta.pdb)}

    P=np.full(Y.shape,np.nan,np.float32); used=np.zeros(len(Y),bool)
    for fold in proto:
        te=np.array([pdb2i[p] for p in fold if p in pdb2i])
        if len(te)<3: continue
        tr=np.setdiff1d(np.arange(len(Y)),te); used[te]=True
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
    print(f"  TEIM-Res  macro-AUPRC={np.nanmean(mm['auprc']):.4f}  macro-AUC={np.nanmean(mm['auc']):.4f}")
    np.savez("results/tables/teim_res_pairwise.npz",
             pdb=meta.pdb.values[idx], clu70=meta.clu70.values[idx],
             auprc=mm["auprc"], auc=mm["auc"], pred=P[idx])
    print("wrote results/tables/teim_res_pairwise.npz")

if __name__=="__main__":
    main()
