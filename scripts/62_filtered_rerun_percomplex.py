#!/usr/bin/env python3
"""Recover per-complex average precision under the SIMILARITY-FILTERED sweep regime, at the three
thresholds that have a size-matched control.

Authorised in results/review/FILTERED_RERUN_AUTHORISATION.md, committed before this ran.
The sweep (scripts/25b_similarity_curve_p82v1.py) stores only aggregates, so the per-complex
values needed for the paired contrast against the size-matched control do not exist anywhere.

The configuration is the sweep's own, copied from that script and not chosen here: the fixed
five-fold partition default_rng(0).permutation(174), the filter Sim[t, te].max() < thr, PairNet,
400 epochs, AdamW lr 1e-3 wd 0.01, BCEWithLogits with the subset's positive-class weight, grad
clip 1.0, seeds 42/43/44 averaged. No tuning, no reselection, no new seeds.

VALIDATION GATE: the reproduced aggregates must match results/tables/similarity_curve_p84v2.csv
to 4 dp. On failure the per-complex values are NOT written for use and the script exits non-zero.

The frozen sweep artifacts are read-only; output goes to results/review/filtered_rerun/.

usage: python scripts/62_filtered_rerun_percomplex.py
"""
import os
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
import sys, time
from pathlib import Path
import numpy as np, pandas as pd, torch, torch.nn as nn
from Levenshtein import ratio

sys.path.insert(0, ".")
from src.data import p84v2 as P84
P84.install_grip_shim()
from src.models.grip import encode_seq, SEG_ORDER, SEG_LENS   # noqa: F401 (used by exec'd source)
from src.eval.metrics import per_complex

_src = open("scripts/19_train_pairmap.py").read().split("def main()")[0]
_src = _src.replace("import argparse, json, sys", "").replace("class PairNet(", "class PairNetCls(")
exec(_src)

# All seven sweep thresholds. 0.95/0.85/0.70 were recovered first for the size-matched control
# (FILTERED_RERUN_AUTHORISATION.md); 1.01/0.90/0.80/0.75 were added to report one authoritative
# uncertainty calculation across the whole sweep (SWEEP_UNCERTAINTY_RECONCILIATION.md).
THRESHOLDS = [1.01, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70]
SEEDS = [42, 43, 44]                   # the sweep's own
EPOCHS = 400
OUT = Path("results/review/filtered_rerun"); OUT.mkdir(parents=True, exist_ok=True)
FROZEN = "results/tables/similarity_curve_p84v2.csv"


def main():
    t0 = time.time()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {dev}", flush=True)
    d = P84.load(); Y, S, pdbs = d["Y"], d["S"], d["pdb"]
    md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
    meta = md.loc[pdbs]
    T, M = P84.T_for(pdbs), P84.M_for(pdbs)
    seqs = meta.tcr_seq.tolist(); n = len(seqs)
    Sim = np.eye(n)
    for i in range(n):
        for j in range(i+1, n):
            Sim[i, j] = Sim[j, i] = ratio(seqs[i], seqs[j])
    folds = np.array_split(np.random.default_rng(0).permutation(n), 5)
    clu = meta.clu70.values
    pdbs = np.array(pdbs)

    frozen = pd.read_csv(FROZEN).set_index("threshold")
    rows, agg = [], []
    done = pd.read_csv(OUT/"filtered_per_complex.tsv", sep="\t") if (OUT/"filtered_per_complex.tsv").exists() else None
    for thr in THRESHOLDS:
        if (OUT/f"filtered_predictions_thr{thr}.npz").exists() and done is not None and (done.threshold == thr).any():
            prev = done[done.threshold == thr]
            model, freq = float(prev.model_auprc.mean()), float(prev.freq_auprc.mean())
            fz = frozen.loc[thr]
            ok = (round(model,4)==round(float(fz.model),4) and round(freq,4)==round(float(fz.freq),4)
                  and round(model-freq,4)==round(float(fz["diff"]),4))
            agg.append(dict(threshold=thr, model=round(model,4), freq=round(freq,4),
                            diff=round(model-freq,4), frozen_model=fz.model, frozen_freq=fz.freq,
                            frozen_diff=fz["diff"], matches_frozen=ok))
            rows.extend(prev.to_dict("records"))
            print(f"  thr {thr}: reusing already-recovered arrays  {'MATCH' if ok else 'MISMATCH'}", flush=True)
            continue
        pm = np.full(Y.shape, np.nan, np.float32); pf = np.full(Y.shape, np.nan, np.float32)
        kept = {}
        for fi, te in enumerate(folds):
            tr = np.setdiff1d(np.arange(n), te)
            tr2 = tr[np.array([Sim[t, te].max() < thr for t in tr])]
            assert len(tr2) >= 30 and not np.isin(tr2, te).any()
            kept[fi] = tr2
            f = (Y[tr2]*S[tr2]).sum(0)/np.maximum(S[tr2].sum(0), 1)
            assert not np.isnan(f).any()
            pf[te] = f[None]
            ps = []
            for sd in SEEDS:
                torch.manual_seed(sd)
                net = PairNetCls(mhc_len=M.shape[1]).to(dev)
                Xa, Ma, Ya, Sa = (torch.tensor(z[tr2]).to(dev) for z in (T, M, Y, S))
                pw = ((Ya == 0) & (Sa > 0)).sum()/torch.clamp(((Ya == 1) & (Sa > 0)).sum(), min=1)
                opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=0.01)
                lf = nn.BCEWithLogitsLoss(reduction="none", pos_weight=pw)
                for _ in range(EPOCHS):
                    net.train(); opt.zero_grad()
                    l = (lf(net(Xa, Ma), Ya)*Sa).sum()/Sa.sum(); l.backward()
                    nn.utils.clip_grad_norm_(net.parameters(), 1.0); opt.step()
                net.eval()
                with torch.no_grad():
                    ps.append(torch.sigmoid(net(torch.tensor(T[te]).to(dev),
                                                torch.tensor(M[te]).to(dev))).cpu().numpy())
            pm[te] = np.mean(ps, 0)
        mm = per_complex(Y, pm, S); ff = per_complex(Y, pf, S)
        model, freq = float(np.nanmean(mm["auprc"])), float(np.nanmean(ff["auprc"]))
        fz = frozen.loc[thr]
        ok = (round(model, 4) == round(float(fz.model), 4)
              and round(freq, 4) == round(float(fz.freq), 4)
              and round(model-freq, 4) == round(float(fz["diff"]), 4))
        agg.append(dict(threshold=thr, model=round(model, 4), freq=round(freq, 4),
                        diff=round(model-freq, 4), frozen_model=fz.model, frozen_freq=fz.freq,
                        frozen_diff=fz["diff"], matches_frozen=ok))
        print(f"  thr {thr}: model={model:.4f} (frozen {fz.model}) freq={freq:.4f} "
              f"(frozen {fz.freq}) diff={model-freq:+.4f} (frozen {fz['diff']:+.4f})  "
              f"{'MATCH' if ok else 'MISMATCH'}", flush=True)
        for i in range(n):
            rows.append(dict(threshold=thr, pdb=pdbs[i], clu70=clu[i],
                             model_auprc=mm["auprc"][i], freq_auprc=ff["auprc"][i],
                             n_pos=mm["n_pos"][i]))
        np.savez_compressed(OUT/f"filtered_predictions_thr{thr}.npz",
                            pred_model=pm, pred_freq=pf, Y=Y, S=S, pdb=pdbs,
                            **{f"kept_fold{k}": pdbs[v] for k, v in kept.items()})

    a = pd.DataFrame(agg); a.to_csv(OUT/"validation_against_frozen.tsv", sep="\t", index=False)
    print("\n=== VALIDATION GATE ===" ); print(a.to_string(index=False))
    if not a.matches_frozen.all():
        print("\nGATE FAILED: reproduced aggregates differ from the frozen sweep. "
              "Per-complex values are NOT written for use; the paired contrast is not computable.",
              flush=True)
        (OUT/"GATE_FAILED").write_text("reproduced aggregates differ from similarity_curve_p84v2.csv\n")
        sys.exit(1)
    pd.DataFrame(rows).to_csv(OUT/"filtered_per_complex.tsv", sep="\t", index=False)
    print(f"\nGATE PASSED. wrote {OUT/'filtered_per_complex.tsv'}  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
