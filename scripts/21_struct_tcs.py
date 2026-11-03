#!/usr/bin/env python3
"""
Struct-TCS / Seq-TCS baseline, reimplemented from Gao et al. (JCIM 2025).

Architecture per the published description:
  MHC branch : Linear -> LSTM -> stacked self-attention (pepSeq blocks)
  TCR branch : Linear -> LSTM -> stacked gated multi-head GAT over residue graph
  Fusion     : cross-attention (TCR query, MHC key/value) + FFN
  Output     : per-residue sigmoid

Reimplemented because the released code (a) omits opt.step() in train_epoch and
(b) calls Multihead_GAT with 6 positional args against a 5-arg signature.
We additionally use a held-out validation fold for early stopping rather than
the test set. Seq-TCS = identical, with the GAT replaced by self-attention.
"""
import argparse, json, sys
import numpy as np, pandas as pd, torch, torch.nn as nn
from sklearn.metrics import roc_auc_score, average_precision_score
sys.path.insert(0,'.')
from src.models.grip import encode_seq, VOCAB, SEG_ORDER, SEG_LENS, TCR_LEN

class GatedGAT(nn.Module):
    """Multi-head graph attention with per-head gating (their Multihead_GaAtN)."""
    def __init__(self, d, heads):
        super().__init__()
        self.h = heads
        self.W = nn.ModuleList([nn.Linear(d, d, bias=False) for _ in range(heads)])
        self.a = nn.ModuleList([nn.Linear(2*d, 1, bias=False) for _ in range(heads)])
        self.gate = nn.Linear(d, heads)
        self.w0 = nn.Linear(d*heads, d, bias=False)
    def forward(self, h, adj):
        B, N, D = h.shape
        g = torch.sigmoid(self.gate(h))                       # B,N,heads
        outs = []
        for k in range(self.h):
            wh = self.W[k](h)
            e = self.a[k](torch.cat([wh.unsqueeze(2).expand(B,N,N,D),
                                     wh.unsqueeze(1).expand(B,N,N,D)], -1)).squeeze(-1)
            e = torch.nn.functional.leaky_relu(e, 0.2).masked_fill(adj == 0, -1e9)
            outs.append(g[:,:,k:k+1] * torch.relu(e.softmax(-1) @ wh))
        return self.w0(torch.cat(outs, -1))

class AttnBlock(nn.Module):
    def __init__(self, d, heads, dropout=0.1):
        super().__init__()
        self.att = nn.MultiheadAttention(d, heads, batch_first=True, dropout=dropout)
        self.n1, self.n2 = nn.LayerNorm(d), nn.LayerNorm(d)
        self.ff = nn.Sequential(nn.Linear(d, d*2), nn.ReLU(), nn.Linear(d*2, d))
    def forward(self, x):
        a,_ = self.att(x,x,x)
        x = self.n1(x+a)
        return self.n2(self.ff(x)+x)

class GATBlock(nn.Module):
    def __init__(self, d, heads):
        super().__init__()
        self.gat = GatedGAT(d, heads)
        self.n1, self.n2 = nn.LayerNorm(d), nn.LayerNorm(d)
        self.ff = nn.Sequential(nn.Linear(d, d*2), nn.ReLU(), nn.Linear(d*2, d))
    def forward(self, x, adj):
        t = self.n1(self.gat(x, adj) + x)
        return self.n2(self.ff(t) + t)

class OursPerRes(nn.Module):
    """Our architecture, per-residue head, for matched comparison."""
    def __init__(self, d=96, heads=6, layers=2, dropout=0.25, mhc_len=28):
        super().__init__()
        self.emb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.seg = nn.Embedding(len(SEG_ORDER), d)
        self.pos = nn.Parameter(torch.randn(1, TCR_LEN, d)*.02)
        self.memb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.mpos = nn.Parameter(torch.randn(1, mhc_len, d)*.02)
        self.enc = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(d, heads, d*2, dropout, batch_first=True, norm_first=True), layers)
        self.menc = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(d, heads, d*2, dropout, batch_first=True, norm_first=True), 1)
        self.cross = nn.MultiheadAttention(d, heads, dropout=dropout, batch_first=True)
        self.norm = nn.LayerNorm(d)
        self.head = nn.Sequential(nn.Linear(d, d), nn.GELU(), nn.Dropout(dropout), nn.Linear(d,1))
        sid,o = torch.zeros(TCR_LEN,dtype=torch.long),0
        for i,ss in enumerate(SEG_ORDER):
            sid[o:o+SEG_LENS[ss]]=i; o+=SEG_LENS[ss]
        self.register_buffer("sid", sid)
    def forward(self, tcr, mhc, adj=None):
        x = self.enc(self.emb(tcr) + self.seg(self.sid)[None] + self.pos)
        m = self.menc(self.memb(mhc) + self.mpos)
        a,_ = self.cross(x, m, m)
        return self.head(self.norm(x + a)).squeeze(-1)

class TCS(nn.Module):
    def __init__(self, d=64, pep_heads=4, prot_heads=4, peplayer=2, protlayer=2,
                 caheads=4, mode="struct", mhc_len=28):
        super().__init__()
        self.mode = mode
        self.emb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.memb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.tlstm = nn.LSTM(d, d, 2, batch_first=True)
        self.mlstm = nn.LSTM(d, d, 2, batch_first=True)
        self.mhc_blocks = nn.ModuleList([AttnBlock(d, pep_heads) for _ in range(peplayer)])
        if mode == "struct":
            self.tcr_blocks = nn.ModuleList([GATBlock(d, prot_heads) for _ in range(protlayer)])
        else:
            self.tcr_blocks = nn.ModuleList([AttnBlock(d, prot_heads) for _ in range(protlayer)])
        self.ca = nn.MultiheadAttention(d, caheads, batch_first=True)
        self.ffn = nn.Sequential(nn.Linear(d, d*2), nn.ReLU(), nn.Linear(d*2, d), nn.LayerNorm(d))
        self.final = nn.Linear(d, 1)
    def forward(self, tcr, mhc, adj):
        m = self.mlstm(self.memb(mhc))[0]
        for b in self.mhc_blocks: m = b(m)
        x = self.tlstm(self.emb(tcr))[0]
        for b in self.tcr_blocks:
            x = b(x, adj) if self.mode == "struct" else b(x)
        inter,_ = self.ca(x, m, m)
        return self.final(self.ffn(inter + x)).squeeze(-1)

def pad_lab(s,L):
    v=np.zeros(L,dtype=np.float32)
    for i,c in enumerate(s[:L]): v[i]=int(c)
    return v

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="struct", choices=["struct","seq","ours"])
    ap.add_argument("--split_by", default="allele", choices=["random","allele","trbv"])
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--epochs", type=int, default=300)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    dev='cuda' if torch.cuda.is_available() else 'cpu'
    Z = np.load("data/processed/pair_maps.npz")
    A = np.load("data/processed/tcr_graphs.npz")["A"]
    meta = pd.read_parquet("data/processed/pair_meta.parquet")
    # per-residue labels: does this CDR residue contact ANY MHC contact position
    Y = (Z["Y"].sum(-1) > 0).astype(np.float32)
    S = (Z["S"].sum(-1) > 0).astype(np.float32)
    mhc_len = Z["Y"].shape[2]
    T = np.stack([np.concatenate([encode_seq(r[k], SEG_LENS[k]) for k in SEG_ORDER]) for _,r in meta.iterrows()])
    M = np.stack([encode_seq(s, mhc_len) for s in meta.mhc_contact])
    Tt,Mt,At = torch.tensor(T), torch.tensor(M), torch.tensor(A)
    Yt,St = torch.tensor(Y), torch.tensor(S)

    groups = meta[a.split_by].fillna("NA").values if a.split_by!="random" else np.arange(len(meta))
    uniq = pd.unique(groups); rng=np.random.default_rng(a.seed); rng.shuffle(uniq)
    folds = np.array_split(uniq, a.folds)
    ys, ps = [], []
    for fi, hold in enumerate(folds):
        te = np.isin(groups,hold); trv = ~te
        if te.sum()<3 or trv.sum()<20: continue
        vi = rng.random(trv.sum()) < 0.15
        tr_idx = np.where(trv)[0][~vi]; va_idx = np.where(trv)[0][vi]
        torch.manual_seed(a.seed)
        net = (OursPerRes(mhc_len=mhc_len) if a.mode=="ours" else TCS(mode=a.mode, mhc_len=mhc_len)).to(dev)
        opt = torch.optim.Adam(net.parameters(), lr=1e-3, weight_decay=1e-4)
        g = lambda i: (Tt[i].to(dev), Mt[i].to(dev), At[i].to(dev), Yt[i].to(dev), St[i].to(dev))
        Xa,Ma,Aa,Ya,Sa = g(tr_idx); Xv,Mv,Av,Yv,Sv = g(va_idx); Xe,Me,Ae,Ye,Se = g(np.where(te)[0])
        pw = ((Ya==0)&(Sa>0)).sum()/torch.clamp(((Ya==1)&(Sa>0)).sum(),min=1)
        lf = nn.BCEWithLogitsLoss(reduction='none', pos_weight=pw)
        best,bad,bs = -1,0,None
        for ep in range(a.epochs):
            net.train(); opt.zero_grad()
            l = (lf(net(Xa,Ma,Aa), Ya)*Sa).sum()/Sa.sum()
            l.backward(); nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
            if ep % 5 == 0:
                net.eval()
                with torch.no_grad(): pv = torch.sigmoid(net(Xv,Mv,Av))
                k = Sv.bool()
                sc = average_precision_score(Yv[k].cpu(), pv[k].cpu())
                if sc > best: best,bad,bs = sc,0,{kk:v.clone() for kk,v in net.state_dict().items()}
                else:
                    bad += 1
                    if bad > 10: break
        net.load_state_dict(bs); net.eval()
        with torch.no_grad(): pe = torch.sigmoid(net(Xe,Me,Ae))
        k = Se.bool()
        ys.append(Ye[k].cpu().numpy()); ps.append(pe[k].cpu().numpy())
    y,p = np.concatenate(ys), np.concatenate(ps)
    auc,apr = roc_auc_score(y,p), average_precision_score(y,p)
    name = {"struct":"Struct-TCS","seq":"Seq-TCS","ours":"Ours"}[a.mode]
    print(f"  {name:11s} [{a.split_by}] seed{a.seed}  AUC={auc:.4f}  AUPRC={apr:.4f}")
    open("results/tables/tcs_runs.jsonl","a").write(json.dumps(dict(
        model=name, split_by=a.split_by, seed=a.seed,
        auc=round(float(auc),4), auprc=round(float(apr),4)))+"\n")

if __name__=="__main__":
    main()
