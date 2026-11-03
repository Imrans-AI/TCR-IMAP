#!/usr/bin/env python3
"""p82v1: pairwise task, all five frozen protocols, GRIP vs. within-fold frequency prior.
Direct copy of scripts/24_run_all.py logic, pointed at the 82-position data, with
per-protocol resume-safe caching.
"""
import argparse, json, sys, time
import numpy as np, pandas as pd, torch, torch.nn as nn
from pathlib import Path
sys.path.insert(0,'.')
from src.models.grip import encode_seq, SEG_ORDER, SEG_LENS, TCR_LEN
from src.eval.metrics import per_complex, cluster_bootstrap_diff
from src.eval.resume import resumable_unit, refuse_overwrite
_src = open("scripts/19_train_pairmap.py").read().split("def main()")[0]
_src = _src.replace("import argparse, json, sys", "").replace("class PairNet(", "class PairNetCls(")
exec(_src)

CACHE_DIR = "results/tables/_partial/p82v1/24b_run_all"
FINAL_CSV = "results/tables/protocol_main_p82v1.csv"

def load():
    Z = np.load("data/processed/pair_maps_p82v1.npz", allow_pickle=True)
    meta = pd.read_parquet("data/processed/pair_meta_p82v1.parquet")
    assert list(Z["pdb"]) == list(meta.pdb)
    md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
    meta = meta.join(md[["clu70","clu80","clu90","dup_family"]], on="pdb")
    T = np.stack([np.concatenate([encode_seq(r[k], SEG_LENS[k]) for k in SEG_ORDER]) for _,r in meta.iterrows()])
    M = np.stack([encode_seq(s, Z["Y"].shape[2]) for s in meta.mhc_contact])
    return Z["Y"].astype(np.float32), Z["S"].astype(np.float32), T, M, meta

def fit_predict(Y,S,T,M,tr,te,seed,epochs=400,dev='cuda'):
    torch.manual_seed(seed)
    net = PairNetCls(mhc_len=M.shape[1]).to(dev)
    Xa,Ma,Ya,Sa = (torch.tensor(z[tr]).to(dev) for z in (T,M,Y,S))
    pw = ((Ya==0)&(Sa>0)).sum()/torch.clamp(((Ya==1)&(Sa>0)).sum(),min=1)
    opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=0.01)
    lf = nn.BCEWithLogitsLoss(reduction='none', pos_weight=pw)
    for ep in range(epochs):
        net.train(); opt.zero_grad()
        l=(lf(net(Xa,Ma),Ya)*Sa).sum()/Sa.sum(); l.backward()
        nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
    net.eval()
    with torch.no_grad():
        return torch.sigmoid(net(torch.tensor(T[te]).to(dev), torch.tensor(M[te]).to(dev))).cpu().numpy()

def run_protocol(pname, folds, Y, S, T, M, meta, dev, seeds):
    preds = {k: np.full(Y.shape, np.nan, np.float32) for k in ("model","pair_freq")}
    covered = np.zeros(len(Y), bool)
    pdb2i = {p:i for i,p in enumerate(meta.pdb)}
    for fold in folds:
        te = np.array([pdb2i[p] for p in fold if p in pdb2i])
        if len(te)<3: continue
        tr = np.setdiff1d(np.arange(len(Y)), te)
        assert set(tr.tolist()) & set(te.tolist()) == set()
        covered[te]=True
        f = (Y[tr]*S[tr]).sum(0)/np.maximum(S[tr].sum(0),1)
        assert not np.isnan(f).any()
        preds["pair_freq"][te] = f[None]
        ps = [fit_predict(Y,S,T,M,tr,te,s,dev=dev) for s in seeds]
        preds["model"][te] = np.mean(ps,0)
    return preds, covered

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--seeds",nargs="+",type=int,default=[42,43,44]); a=ap.parse_args()
    dev='cuda' if torch.cuda.is_available() else 'cpu'
    print(f"device: {dev}")
    refuse_overwrite(FINAL_CSV)
    Y,S,T,M,meta = load()
    proto = json.load(open("data/splits_v2/protocols.json"))
    rows=[]
    for pname, folds in proto.items():
        t0 = time.time()
        def compute():
            preds, covered = run_protocol(pname, folds, Y, S, T, M, meta, dev, tuple(a.seeds))
            return dict(pdb=meta.pdb.values, clu70=meta.clu70.values, covered=covered,
                        Y=Y.astype(np.int8), S=S.astype(np.int8),
                        pred_model=preds["model"], pred_freq=preds["pair_freq"])
        def validate(d):
            # pred_model/pred_freq are legitimately NaN for complexes not covered by any
            # usable fold of this protocol (e.g. the excluded n=2 allele fold) -- only the
            # covered subset must be finite.
            cov = d["covered"].astype(bool)
            if not np.isfinite(d["pred_model"][cov]).all(): return False
            if not np.isfinite(d["pred_freq"][cov]).all(): return False
            return bool(cov.any())
        result, resumed = resumable_unit(CACHE_DIR, f"protocol_{pname}", compute, validate_fn=validate)
        dt = time.time()-t0
        idx = np.where(result["covered"])[0]
        clusters = result["clu70"][idx]
        mm = {k: per_complex(result["Y"][idx], result[f"pred_{k}"][idx], result["S"][idx]) for k in ("model","freq")}
        d,lo,hi,p = cluster_bootstrap_diff(mm["model"]["auprc"], mm["freq"]["auprc"], clusters)
        line = dict(protocol=pname, n=len(idx),
            model_auprc=round(float(np.nanmean(mm["model"]["auprc"])),4),
            freq_auprc=round(float(np.nanmean(mm["freq"]["auprc"])),4),
            model_p_at_L=round(float(np.nanmean(mm["model"]["p_at_L"])),4),
            freq_p_at_L=round(float(np.nanmean(mm["freq"]["p_at_L"])),4),
            diff=round(d,4), ci_lo=round(lo,4), ci_hi=round(hi,4), pval=round(p,4))
        rows.append(line)
        print(f"[{pname}] n={len(idx)} model={line['model_auprc']:.4f} freq={line['freq_auprc']:.4f} "
              f"diff={d:+.4f} [{lo:+.4f},{hi:+.4f}] p={p:.4f} (resumed={resumed}, {dt:.1f}s)")
        pm_path = f"results/tables/pairmaps_{pname}_p82v1.npz"
        if Path(pm_path).exists():
            print(f"  {pm_path} already finalized, skipping rewrite (resume)")
        else:
            pm_tmp = pm_path.replace(".npz", ".tmp.npz")
            np.savez_compressed(pm_tmp, pdb=result["pdb"], clu70=result["clu70"], covered=result["covered"],
                                Y=result["Y"], S=result["S"], pred_model=result["pred_model"], pred_freq=result["pred_freq"])
            Path(pm_tmp).replace(pm_path)
            print(f"  wrote {pm_path}")
    df = pd.DataFrame(rows)
    tmp = FINAL_CSV + ".tmp"
    df.to_csv(tmp, index=False)
    Path(tmp).replace(FINAL_CSV)
    print(f"\nwrote {FINAL_CSV}")
    print(df.to_string(index=False))

if __name__=="__main__":
    main()
