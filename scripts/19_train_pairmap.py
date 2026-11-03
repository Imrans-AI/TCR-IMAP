#!/usr/bin/env python3
"""Pairwise CDR-MHC contact map prediction + matched baselines."""
import argparse, json, sys
import numpy as np, pandas as pd, torch, torch.nn as nn
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.linear_model import LogisticRegression
sys.path.insert(0,'.')
from src.models.grip import encode_seq, VOCAB, SEG_ORDER, SEG_LENS, TCR_LEN

class PairNet(nn.Module):
    """Outputs a full 84x28 interaction map from a bilinear head over paired embeddings."""
    def __init__(self, d=128, heads=8, layers=2, dropout=0.2, mhc_len=28):
        super().__init__()
        self.emb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.seg = nn.Embedding(len(SEG_ORDER), d)
        self.pos = nn.Parameter(torch.randn(1, TCR_LEN, d)*.02)
        self.memb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.mpos = nn.Parameter(torch.randn(1, mhc_len, d)*.02)
        self.tenc = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(d, heads, d*2, dropout, batch_first=True, norm_first=True), layers)
        self.menc = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(d, heads, d*2, dropout, batch_first=True, norm_first=True), 1)
        self.cross = nn.MultiheadAttention(d, heads, dropout=dropout, batch_first=True)
        self.norm = nn.LayerNorm(d)
        self.wt = nn.Linear(d, d, bias=False)      # bilinear interaction
        self.wm = nn.Linear(d, d, bias=False)
        self.pair_pos = nn.Parameter(torch.zeros(TCR_LEN, mhc_len))   # learned positional prior
        self.scale = nn.Parameter(torch.tensor(1.0))
        self.bias = nn.Parameter(torch.tensor(-4.0))
        sid,o = torch.zeros(TCR_LEN,dtype=torch.long),0
        for i,s in enumerate(SEG_ORDER):
            sid[o:o+SEG_LENS[s]]=i; o+=SEG_LENS[s]
        self.register_buffer("sid", sid)
    def forward(self, t, m):
        x = self.tenc(self.emb(t) + self.seg(self.sid)[None] + self.pos)
        mm = self.menc(self.memb(m) + self.mpos)
        a,_ = self.cross(x, mm, mm)
        x = self.norm(x + a)
        logits = torch.einsum('bid,bjd->bij', self.wt(x), self.wm(mm)) / x.shape[-1]**0.5
        return logits*self.scale + self.pair_pos[None] + self.bias

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split_by", default="allele", choices=["random","allele","trbv"])
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--epochs", type=int, default=400)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    dev='cuda' if torch.cuda.is_available() else 'cpu'
    Z = np.load("data/processed/pair_maps.npz"); Y, S = Z["Y"], Z["S"]
    meta = pd.read_parquet("data/processed/pair_meta.parquet")
    mhc_len = Y.shape[2]
    T = np.stack([np.concatenate([encode_seq(r[k], SEG_LENS[k]) for k in SEG_ORDER]) for _,r in meta.iterrows()])
    M = np.stack([encode_seq(s, mhc_len) for s in meta.mhc_contact])
    Tt,Mt = torch.tensor(T), torch.tensor(M)
    Yt,St = torch.tensor(Y,dtype=torch.float32), torch.tensor(S,dtype=torch.float32)

    groups = meta[a.split_by].fillna("NA").values if a.split_by!="random" else np.arange(len(meta))
    uniq = pd.unique(groups); rng=np.random.default_rng(a.seed); rng.shuffle(uniq)
    folds = np.array_split(uniq, a.folds)

    coll = {k:{"y":[],"p":[]} for k in ("model","pair_freq","logreg")}
    for hold in folds:
        te = np.isin(groups,hold); tr = ~te
        if te.sum()<3 or tr.sum()<20: continue
        # baseline 1: empirical pair frequency
        f = (Y[tr]*S[tr]).sum(0)/np.maximum(S[tr].sum(0),1)
        # baseline 2: logreg on (tcr residue, mhc residue, i, j)
        def flat(idx):
            n=idx.sum()
            ti = np.repeat(T[idx][:,:,None], mhc_len, axis=2).reshape(-1)
            mi = np.repeat(M[idx][:,None,:], TCR_LEN, axis=1).reshape(-1)
            ii = np.tile(np.arange(TCR_LEN)[:,None], (1,mhc_len)); ii=np.tile(ii,(n,1,1)).reshape(-1)
            jj = np.tile(np.arange(mhc_len)[None,:], (TCR_LEN,1)); jj=np.tile(jj,(n,1,1)).reshape(-1)
            X = np.concatenate([np.eye(VOCAB)[ti], np.eye(VOCAB)[mi],
                                np.eye(TCR_LEN)[ii], np.eye(mhc_len)[jj]], 1).astype(np.float32)
            return X, Y[idx].reshape(-1), S[idx].reshape(-1)
        Xtr,ytr,str_ = flat(tr); Xte,yte,ste = flat(te)
        sub = np.where(str_>0)[0]
        sub = rng.choice(sub, min(200000,len(sub)), replace=False)
        lr = LogisticRegression(max_iter=400, C=0.5, solver='liblinear').fit(Xtr[sub], ytr[sub])
        k = ste>0
        coll["pair_freq"]["p"].append(np.tile(f,(te.sum(),1)).reshape(-1)[k]); coll["pair_freq"]["y"].append(yte[k])
        coll["logreg"]["p"].append(lr.predict_proba(Xte[k])[:,1]); coll["logreg"]["y"].append(yte[k])
        # model
        torch.manual_seed(a.seed)
        net = PairNet(mhc_len=mhc_len).to(dev)
        Xa,Ma,Ya,Sa = Tt[tr].to(dev), Mt[tr].to(dev), Yt[tr].to(dev), St[tr].to(dev)
        pw = ((Ya==0)&(Sa>0)).sum()/torch.clamp(((Ya==1)&(Sa>0)).sum(),min=1)
        opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=0.01)
        lf = nn.BCEWithLogitsLoss(reduction='none', pos_weight=pw)
        for ep in range(a.epochs):
            net.train(); opt.zero_grad()
            l = (lf(net(Xa,Ma), Ya)*Sa).sum()/Sa.sum()
            l.backward(); nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
        net.eval()
        with torch.no_grad():
            p = torch.sigmoid(net(Tt[te].to(dev), Mt[te].to(dev))).cpu().numpy()
        coll["model"]["p"].append(p.reshape(-1)[k]); coll["model"]["y"].append(yte[k])

    print(f"\n[{a.split_by} split, seed {a.seed}]")
    out={}
    for m,d in coll.items():
        y=np.concatenate(d["y"]); p=np.concatenate(d["p"])
        auc,apr = roc_auc_score(y,p), average_precision_score(y,p)
        out[m]=(auc,apr); print(f"  {m:10s} AUC={auc:.4f}  AUPRC={apr:.4f}")
    print(f"  model - logreg : {out['model'][0]-out['logreg'][0]:+.4f} AUC, {out['model'][1]-out['logreg'][1]:+.4f} AUPRC")
    open("results/tables/pairmap_runs.jsonl","a").write(json.dumps(dict(
        split_by=a.split_by, seed=a.seed,
        **{f"{m}_{k}":round(float(v),4) for m,(x,z) in out.items() for k,v in (("auc",x),("auprc",z))}))+"\n")

if __name__=="__main__":
    main()
