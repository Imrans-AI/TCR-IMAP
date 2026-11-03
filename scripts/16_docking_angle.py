#!/usr/bin/env python3
"""Task 2: docking-angle regression from CDR + MHC contact sequence."""
import argparse, json, sys
import numpy as np, pandas as pd, torch, torch.nn as nn
from scipy.stats import pearsonr, spearmanr
sys.path.insert(0,'.')
from src.models.grip import encode_seq, VOCAB, SEG_ORDER, SEG_LENS, TCR_LEN

class AngleNet(nn.Module):
    def __init__(self, d=96, heads=6, dropout=0.3):
        super().__init__()
        self.emb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.seg = nn.Embedding(len(SEG_ORDER), d)
        self.pos = nn.Parameter(torch.randn(1, TCR_LEN, d)*.02)
        self.memb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.mpos = nn.Parameter(torch.randn(1, 28, d)*.02)
        self.enc = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(d, heads, d*2, dropout, batch_first=True, norm_first=True), 2)
        self.cross = nn.MultiheadAttention(d, heads, dropout=dropout, batch_first=True)
        self.norm = nn.LayerNorm(d)
        self.head = nn.Sequential(nn.Linear(d*2, d), nn.GELU(), nn.Dropout(dropout), nn.Linear(d,1))
        sid,o = torch.zeros(TCR_LEN,dtype=torch.long),0
        for i,s in enumerate(SEG_ORDER):
            sid[o:o+SEG_LENS[s]]=i; o+=SEG_LENS[s]
        self.register_buffer("sid", sid)
    def forward(self, t, m):
        x = self.enc(self.emb(t) + self.seg(self.sid)[None] + self.pos)
        mm = self.memb(m) + self.mpos
        a,_ = self.cross(x, mm, mm)
        x = self.norm(x + a)
        return self.head(torch.cat([x.mean(1), mm.mean(1)], -1)).squeeze(-1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split_by", default="trbv", choices=["random","allele","trbv"])
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--epochs", type=int, default=400)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    df = pd.read_parquet("data/processed/contact_sites.parquet").reset_index(drop=True)
    df["ang"] = pd.to_numeric(df.docking_angle, errors="coerce")
    df = df[df.ang.notna()].reset_index(drop=True)

    T = np.stack([np.concatenate([encode_seq(r[k], SEG_LENS[k]) for k in SEG_ORDER]) for _,r in df.iterrows()])
    M = np.stack([encode_seq(s, 28) for s in df.mhc_contact])
    y = df.ang.values.astype(np.float32)
    T,M = torch.tensor(T), torch.tensor(M)
    Y = torch.tensor(y)

    groups = df[a.split_by].fillna("NA").values if a.split_by!="random" else np.arange(len(df))
    uniq = pd.unique(groups); rng=np.random.default_rng(a.seed); rng.shuffle(uniq)
    folds = np.array_split(uniq, a.folds)
    preds = np.full(len(df), np.nan)
    for hold in folds:
        te = np.isin(groups, hold); tr = ~te
        if te.sum()<3 or tr.sum()<20: continue
        mu, sd = y[tr].mean(), y[tr].std()
        Xtr,Mtr = T[tr].to(dev), M[tr].to(dev)
        ytr = torch.tensor((y[tr]-mu)/sd).to(dev)
        torch.manual_seed(a.seed)
        net = AngleNet().to(dev)
        opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=0.05)
        lf = nn.SmoothL1Loss()
        for ep in range(a.epochs):
            net.train(); opt.zero_grad()
            l = lf(net(Xtr,Mtr), ytr); l.backward()
            nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
        net.eval()
        with torch.no_grad():
            p = net(T[te].to(dev), M[te].to(dev)).cpu().numpy()*sd + mu
        preds[te] = p
    ok = ~np.isnan(preds)
    mae = np.abs(preds[ok]-y[ok]).mean()
    base = np.abs(y[ok]-np.median(y[ok])).mean()
    r_p = pearsonr(preds[ok], y[ok]); r_s = spearmanr(preds[ok], y[ok])
    print(f"  [{a.split_by} split, seed {a.seed}] n={ok.sum()}")
    print(f"    MAE={mae:.2f}deg  (median baseline {base:.2f})  R2={1-((preds[ok]-y[ok])**2).sum()/((y[ok]-y[ok].mean())**2).sum():.3f}")
    print(f"    Pearson r={r_p[0]:.3f} (p={r_p[1]:.2g})   Spearman={r_s[0]:.3f}")
    open("results/tables/angle_runs.jsonl","a").write(json.dumps(dict(
        split_by=a.split_by, seed=a.seed, mae=round(float(mae),3),
        baseline_mae=round(float(base),3), pearson=round(float(r_p[0]),3),
        pval=float(r_p[1]), n=int(ok.sum())))+"\n")

if __name__ == "__main__":
    main()
