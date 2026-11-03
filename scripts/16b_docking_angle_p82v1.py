#!/usr/bin/env python3
"""p82v1: docking-angle regression. Direct copy of scripts/16_docking_angle.py logic,
with two required fixes (per results/review/PRIMARY_82_POSITION_RERUN_PLAN.md section 3):
the two hardcoded `28`s (AngleNet.mpos shape, encode_seq call) are replaced with a real
mhc_len parameter, and the input metadata is data/processed/contact_sites_p82v1.parquet
(82-character mhc_contact, 6uz1 already excluded at build time). R-squared is computed
AND PERSISTED this time, closing the gap flagged in the main audit (section 4.8): the
original 28-position run computed R2 but never saved it to any artifact.

split_by in {random, allele, trbv} matches the ad hoc, non-frozen-protocol split scheme
the original script uses (this analysis was never run under data/splits_v2/protocols.json
-- see the plan's section 3 finding on this point; unchanged here).
"""
import json, sys, time
import numpy as np, pandas as pd, torch, torch.nn as nn
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
sys.path.insert(0,'.')
from src.models.grip import encode_seq, VOCAB, SEG_ORDER, SEG_LENS, TCR_LEN
from src.eval.resume import resumable_unit, refuse_overwrite

CACHE_DIR = "results/tables/_partial/p82v1/16b_docking_angle"
FINAL_JSONL = "results/tables/angle_runs_p82v1.jsonl"

class AngleNet(nn.Module):
    def __init__(self, mhc_len, d=96, heads=6, dropout=0.3):
        super().__init__()
        self.emb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.seg = nn.Embedding(len(SEG_ORDER), d)
        self.pos = nn.Parameter(torch.randn(1, TCR_LEN, d)*.02)
        self.memb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.mpos = nn.Parameter(torch.randn(1, mhc_len, d)*.02)     # was hardcoded 28
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

def run_one(split_by, seed, df, T, M, y, mhc_len, dev, folds=5, epochs=400):
    groups = df[split_by].fillna("NA").values if split_by!="random" else np.arange(len(df))
    uniq = pd.unique(groups); rng=np.random.default_rng(seed); rng.shuffle(uniq)
    fold_groups = np.array_split(uniq, folds)
    preds = np.full(len(df), np.nan)
    for hold in fold_groups:
        te = np.isin(groups, hold); tr = ~te
        if te.sum()<3 or tr.sum()<20: continue
        assert not (te & tr).any()
        mu, sd = y[tr].mean(), y[tr].std()
        Xtr,Mtr = torch.tensor(T[tr]).to(dev), torch.tensor(M[tr]).to(dev)
        ytr = torch.tensor((y[tr]-mu)/sd).to(dev)
        torch.manual_seed(seed)
        net = AngleNet(mhc_len=mhc_len).to(dev)
        opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=0.05)
        lf = nn.SmoothL1Loss()
        for ep in range(epochs):
            net.train(); opt.zero_grad()
            l = lf(net(Xtr,Mtr), ytr); l.backward()
            nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
        net.eval()
        with torch.no_grad():
            p = net(torch.tensor(T[te]).to(dev), torch.tensor(M[te]).to(dev)).cpu().numpy()*sd + mu
        preds[te] = p
    ok = ~np.isnan(preds)
    mae = np.abs(preds[ok]-y[ok]).mean()
    base = np.abs(y[ok]-np.median(y[ok])).mean()
    r_p = pearsonr(preds[ok], y[ok]); r_s = spearmanr(preds[ok], y[ok])
    r2 = 1 - ((preds[ok]-y[ok])**2).sum() / ((y[ok]-y[ok].mean())**2).sum()
    return dict(mae=float(mae), baseline_mae=float(base), pearson=float(r_p[0]),
                pval=float(r_p[1]), spearman=float(r_s[0]), r2=float(r2), n=int(ok.sum()))

def main():
    t_all = time.time()
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"device: {dev}")
    refuse_overwrite(FINAL_JSONL)
    df = pd.read_parquet("data/processed/contact_sites_p82v1.parquet").reset_index(drop=True)
    assert "6uz1" not in df.pdb.values
    df["ang"] = pd.to_numeric(df.docking_angle, errors="coerce")
    df = df[df.ang.notna()].reset_index(drop=True)
    print(f"complexes with a valid docking angle (6uz1 excluded, primary-set universe): n={len(df)}")

    mhc_len = df.mhc_contact.str.len().iloc[0]
    assert (df.mhc_contact.str.len()==mhc_len).all() and mhc_len==82
    T = np.stack([np.concatenate([encode_seq(r[k], SEG_LENS[k]) for k in SEG_ORDER]) for _,r in df.iterrows()])
    M = np.stack([encode_seq(s, mhc_len) for s in df.mhc_contact])
    y = df.ang.values.astype(np.float32)

    rows = []
    for split_by in ("random","allele","trbv"):
        for seed in (42,43,44):
            unit_id = f"{split_by}_seed{seed}"
            t0=time.time()
            def compute():
                r = run_one(split_by, seed, df, T, M, y, mhc_len, dev)
                return {k: np.asarray(v) for k,v in r.items()}
            result, resumed = resumable_unit(CACHE_DIR, unit_id, compute)
            row = {k: (float(v) if v.dtype.kind=='f' else int(v)) for k,v in result.items()}
            row.update(split_by=split_by, seed=seed)
            rows.append(row)
            print(f"  [{split_by} split, seed {seed}] n={row['n']} MAE={row['mae']:.2f} "
                  f"(baseline {row['baseline_mae']:.2f}) r={row['pearson']:.3f} R2={row['r2']:.3f} "
                  f"(resumed={resumed}, {time.time()-t0:.1f}s)")

    tmp = FINAL_JSONL + ".tmp"
    with open(tmp, "w") as f:
        for r in rows:
            f.write(json.dumps(dict(split_by=r["split_by"], seed=r["seed"],
                mae=round(r["mae"],3), baseline_mae=round(r["baseline_mae"],3),
                pearson=round(r["pearson"],3), pval=r["pval"], spearman=round(r["spearman"],3),
                r2=round(r["r2"],4), n=r["n"])) + "\n")
    Path(tmp).replace(FINAL_JSONL)
    print(f"\nwrote {FINAL_JSONL}  (total {time.time()-t_all:.1f}s)")

    dfres = pd.DataFrame(rows)
    print(dfres.groupby("split_by").agg(mae_mean=("mae","mean"), mae_min=("mae","min"), mae_max=("mae","max"),
          r_min=("pearson","min"), r_max=("pearson","max"),
          r2_min=("r2","min"), r2_max=("r2","max")).round(3).to_string())

if __name__=="__main__":
    main()
