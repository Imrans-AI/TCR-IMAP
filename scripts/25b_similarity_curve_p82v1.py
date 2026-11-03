#!/usr/bin/env python3
"""p82v1: identity-threshold sweep (Fig. 3A/B's data source). Direct copy of
scripts/25_similarity_curve.py logic, pointed at the 82-position data. Same fixed
THRESHOLDS list, same deterministic seed=0 internal 5-fold split (independent of
data/splits_v2/protocols.json, exactly matching the original script's own design),
with per-threshold resume-safe caching.
"""
import json, sys, time
import numpy as np, pandas as pd, torch, torch.nn as nn
from pathlib import Path
from Levenshtein import ratio
sys.path.insert(0,'.')
from src.models.grip import encode_seq, SEG_ORDER, SEG_LENS
from src.eval.metrics import per_complex, cluster_bootstrap_diff
from src.eval.resume import resumable_unit, refuse_overwrite

_src = open("scripts/19_train_pairmap.py").read().split("def main()")[0]
_src = _src.replace("import argparse, json, sys","").replace("class PairNet(","class PairNetCls(")
exec(_src)

THRESHOLDS = [1.01, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70]
SEEDS = [42, 43, 44]
CACHE_DIR = "results/tables/_partial/p82v1/25b_similarity_curve"
FINAL_CSV = "results/tables/similarity_curve_p82v1.csv"

def run_threshold(thr, Sim, folds, Y, S, T, M, meta, n, dev):
    preds_m = np.full(Y.shape, np.nan, np.float32)
    preds_f = np.full(Y.shape, np.nan, np.float32)
    used = np.zeros(n, bool)
    for te in folds:
        tr = np.setdiff1d(np.arange(n), te)
        keep = np.array([Sim[t, te].max() < thr for t in tr])
        tr2 = tr[keep]
        if len(tr2) < 30: continue
        assert not np.isin(tr2, te).any()
        used[te] = True
        f = (Y[tr2]*S[tr2]).sum(0)/np.maximum(S[tr2].sum(0),1)
        assert not np.isnan(f).any()
        preds_f[te] = f[None]
        ps = []
        for sd in SEEDS:
            torch.manual_seed(sd)
            net = PairNetCls(mhc_len=M.shape[1]).to(dev)
            Xa,Ma,Ya,Sa = (torch.tensor(z[tr2]).to(dev) for z in (T,M,Y,S))
            pw = ((Ya==0)&(Sa>0)).sum()/torch.clamp(((Ya==1)&(Sa>0)).sum(),min=1)
            opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=0.01)
            lf = nn.BCEWithLogitsLoss(reduction='none', pos_weight=pw)
            for ep in range(400):
                net.train(); opt.zero_grad()
                l=(lf(net(Xa,Ma),Ya)*Sa).sum()/Sa.sum(); l.backward()
                nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
            net.eval()
            with torch.no_grad():
                ps.append(torch.sigmoid(net(torch.tensor(T[te]).to(dev),
                                            torch.tensor(M[te]).to(dev))).cpu().numpy())
        preds_m[te] = np.mean(ps,0)
    idx = np.where(used)[0]
    mm = per_complex(Y[idx], preds_m[idx], S[idx])
    ff = per_complex(Y[idx], preds_f[idx], S[idx])
    d,lo,hi,p = cluster_bootstrap_diff(mm["auprc"], ff["auprc"], meta.clu70.values[idx])
    realised = np.mean([Sim[i, np.setdiff1d(np.arange(n),[i])][Sim[i, np.setdiff1d(np.arange(n),[i])] < thr].max()
                        if (Sim[i, np.setdiff1d(np.arange(n),[i])] < thr).any() else 0 for i in idx])
    return dict(threshold=thr, n=len(idx),
                model=float(np.nanmean(mm["auprc"])), freq=float(np.nanmean(ff["auprc"])),
                diff=float(d), lo=float(lo), hi=float(hi), p=float(p), mean_max_sim=float(realised))

def main():
    t_all = time.time()
    dev='cuda' if torch.cuda.is_available() else 'cpu'
    print(f"device: {dev}")
    refuse_overwrite(FINAL_CSV)
    Z = np.load("data/processed/pair_maps_p82v1.npz", allow_pickle=True)
    Y, S = Z["Y"].astype(np.float32), Z["S"].astype(np.float32)
    meta = pd.read_parquet("data/processed/pair_meta_p82v1.parquet")
    assert list(Z["pdb"]) == list(meta.pdb)
    md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
    meta = meta.join(md[["clu70","tcr_seq"]], on="pdb")
    T = np.stack([np.concatenate([encode_seq(r[k], SEG_LENS[k]) for k in SEG_ORDER]) for _,r in meta.iterrows()])
    M = np.stack([encode_seq(s, Y.shape[2]) for s in meta.mhc_contact])
    seqs = meta.tcr_seq.tolist(); n = len(seqs)

    print("computing pairwise identity matrix ...", flush=True)
    Sim = np.eye(n)
    for i in range(n):
        for j in range(i+1, n):
            Sim[i,j] = Sim[j,i] = ratio(seqs[i], seqs[j])

    rng = np.random.default_rng(0)
    folds = np.array_split(rng.permutation(n), 5)

    rows = []
    for thr in THRESHOLDS:
        t0 = time.time()
        def compute():
            r = run_threshold(thr, Sim, folds, Y, S, T, M, meta, n, dev)
            return {k: np.asarray(v) for k,v in r.items()}
        result, resumed = resumable_unit(CACHE_DIR, f"thr_{thr}", compute)
        row = {k: (v.item() if v.ndim==0 else v) for k,v in result.items()}
        rows.append(row)
        print(f"  thr<{thr:<5} n={int(row['n']):3d} model={row['model']:.4f} freq={row['freq']:.4f} "
              f"diff={row['diff']:+.4f} [{row['lo']:+.4f},{row['hi']:+.4f}] p={row['p']:.4f} "
              f"(resumed={resumed}, {time.time()-t0:.1f}s)")

    df = pd.DataFrame(rows)
    for c in ("threshold","model","freq","diff","lo","hi","p","mean_max_sim"): df[c] = df[c].round(4 if c!="threshold" else 2)
    df["n"] = df["n"].astype(int)
    tmp = FINAL_CSV + ".tmp"
    df.to_csv(tmp, index=False)
    Path(tmp).replace(FINAL_CSV)
    print(f"\nwrote {FINAL_CSV}  (total {time.time()-t_all:.1f}s)")

if __name__=="__main__":
    main()
