#!/usr/bin/env python3
"""Training-size-matched control for the identity-threshold sweep.

Specified in results/review/SIZE_MATCHED_CONTROL_PLAN.md, committed before this file was written.
Nothing here may be changed to chase an outcome; the thresholds, draw seeds, training seeds,
metric and stopping rules are fixed by that document.

The sweep (scripts/25b_similarity_curve_p82v1.py under the p84v2 runner) holds the test set fixed
and removes from training every complex whose similarity to any test complex of the fold reaches
the threshold. That shrinks training from 139-140 complexes per fold to 69-78 at 0.70, so the
sweep cannot separate homology control from training size. This control trains on the SAME NUMBER
of complexes drawn at random from the unrestricted pool, with no similarity filter, on the same
test complexes -- so the only difference from the sweep is which complexes were removed, not how
many.

Held exactly as in the sweep: the fixed five-fold partition default_rng(0).permutation(174),
all 174 test complexes, the p84v2 representation, PairNet from 19_train_pairmap, 400 epochs full
batch, AdamW lr 1e-3 wd 0.01, BCEWithLogits with the subset's positive-class weight, grad clip 1.0,
training seeds 42/43/44 averaged, and the pair-frequency prior fitted on exactly the same subset.

usage: python scripts/60_size_matched_control.py [--check | --run]
  --check  builds inputs, draws every subset, writes the membership and achieved-similarity
           tables, instantiates the model -- no training.
  --run    trains all 225 fits and writes the metric tables.
"""
import os
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
import argparse, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd, torch, torch.nn as nn
from Levenshtein import ratio

sys.path.insert(0, ".")
from src.data import p84v2 as P84
P84.install_grip_shim()                      # PairNet must be built for the 75 x 84 axes
from src.models.grip import encode_seq, SEG_ORDER, SEG_LENS   # noqa: F401  (used by the exec'd source)
from src.eval.metrics import per_complex, cluster_bootstrap_diff

# PairNet exactly as the sweep gets it: the consumer's own exec of 19_train_pairmap's header.
_src = open("scripts/19_train_pairmap.py").read().split("def main()")[0]
_src = _src.replace("import argparse, json, sys", "").replace("class PairNet(", "class PairNetCls(")
exec(_src)

# ---- frozen by results/review/SIZE_MATCHED_CONTROL_PLAN.md section 2 --------------------------
THRESHOLDS = [0.95, 0.85, 0.70]
DRAW_SEEDS = [1001, 1002, 1003, 1004, 1005]
TRAIN_SEEDS = [42, 43, 44]
EPOCHS = 400
N_BOOT, BOOT_SEED = 20000, 42
OUT = Path("results/review/size_matched_control")
SWEEP_CSV = "results/tables/similarity_curve_p84v2.csv"   # read-only input


def load():
    d = P84.load()
    Y, S = d["Y"], d["S"]
    pdbs = d["pdb"]
    md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
    meta = md.loc[pdbs]
    T = P84.T_for(pdbs); M = P84.M_for(pdbs)
    seqs = meta.tcr_seq.tolist(); n = len(seqs)
    Sim = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            Sim[i, j] = Sim[j, i] = ratio(seqs[i], seqs[j])
    # the sweep's own partition -- NOT protocols.json
    folds = np.array_split(np.random.default_rng(0).permutation(n), 5)
    return Y, S, T, M, meta, np.array(pdbs), Sim, folds, n


def retained_counts(Sim, folds, n):
    """k per (threshold, fold): how many complexes the sweep's similarity filter keeps."""
    k = {}
    for thr in THRESHOLDS:
        for fi, te in enumerate(folds):
            tr = np.setdiff1d(np.arange(n), te)
            k[(thr, fi)] = int(np.array([Sim[t, te].max() < thr for t in tr]).sum())
    return k


def draw_subsets(k, folds, n):
    """Size-matched random subsets. One default_rng(seed) per draw, consumed over
    (threshold, fold) in the order fixed by the plan: thresholds in THRESHOLDS order,
    folds 0..4 within each. Declared in advance so the sequence cannot be re-rolled."""
    subs = {}
    for ds in DRAW_SEEDS:
        rng = np.random.default_rng(ds)
        for thr in THRESHOLDS:
            for fi, te in enumerate(folds):
                tr = np.setdiff1d(np.arange(n), te)
                subs[(ds, thr, fi)] = np.sort(rng.choice(tr, k[(thr, fi)], replace=False))
    return subs


def fit_predict(T, M, Y, S, sub, te, dev):
    """The sweep's training loop verbatim, on the sampled subset."""
    ps = []
    for sd in TRAIN_SEEDS:
        torch.manual_seed(sd)
        net = PairNetCls(mhc_len=M.shape[1]).to(dev)
        Xa, Ma, Ya, Sa = (torch.tensor(z[sub]).to(dev) for z in (T, M, Y, S))
        pw = ((Ya == 0) & (Sa > 0)).sum() / torch.clamp(((Ya == 1) & (Sa > 0)).sum(), min=1)
        opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=0.01)
        lf = nn.BCEWithLogitsLoss(reduction="none", pos_weight=pw)
        for _ in range(EPOCHS):
            net.train(); opt.zero_grad()
            l = (lf(net(Xa, Ma), Ya) * Sa).sum() / Sa.sum(); l.backward()
            nn.utils.clip_grad_norm_(net.parameters(), 1.0); opt.step()
        net.eval()
        with torch.no_grad():
            ps.append(torch.sigmoid(net(torch.tensor(T[te]).to(dev),
                                        torch.tensor(M[te]).to(dev))).cpu().numpy())
    return np.mean(ps, 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {dev}", flush=True)

    Y, S, T, M, meta, pdbs, Sim, folds, n = load()
    assert Y.shape[1:] == (75, 84), Y.shape
    k = retained_counts(Sim, folds, n)
    subs = draw_subsets(k, folds, n)

    # membership + achieved similarity, written before any training
    rows = []
    for (ds, thr, fi), sub in subs.items():
        te = folds[fi]
        tr = np.setdiff1d(np.arange(n), te)
        keep = tr[np.array([Sim[t, te].max() < thr for t in tr])]
        rows.append(dict(draw_seed=ds, threshold=thr, fold=fi,
                         n_test=len(te), n_pool=len(tr), n_train=len(sub),
                         n_train_sweep=len(keep),
                         achieved_max_sim=round(float(max(Sim[s, te].max() for s in sub)), 4),
                         achieved_mean_max_sim=round(float(np.mean([Sim[s, te].max() for s in sub])), 4),
                         sweep_max_sim=round(float(max(Sim[s, te].max() for s in keep)), 4),
                         overlap_with_sweep_subset=int(len(np.intersect1d(sub, keep))),
                         train_pdbs=";".join(pdbs[sub])))
    mem = pd.DataFrame(rows).sort_values(["threshold", "draw_seed", "fold"])
    mem.to_csv(OUT/"memberships.tsv", sep="\t", index=False)
    print(f"wrote {OUT/'memberships.tsv'}  ({len(mem)} subsets)", flush=True)
    print(mem.groupby("threshold")[["n_train", "n_train_sweep", "achieved_max_sim",
                                    "sweep_max_sim", "achieved_mean_max_sim"]].mean().round(4).to_string())
    if not a.run:
        net = PairNetCls(mhc_len=M.shape[1])
        print(f"\n--check: PairNet instantiated, {sum(p.numel() for p in net.parameters()):,} parameters")
        print(f"--check: {len(subs)} subsets drawn, {len(subs)*len(TRAIN_SEEDS)} fits would run. No training.")
        return

    clu = meta.clu70.values
    per_rows, draw_rows = [], []
    t_all = time.time()
    for thr in THRESHOLDS:
        store = {}
        for ds in DRAW_SEEDS:
            t0 = time.time()
            pm = np.full(Y.shape, np.nan, np.float32)
            pf = np.full(Y.shape, np.nan, np.float32)
            for fi, te in enumerate(folds):
                sub = subs[(ds, thr, fi)]
                assert not np.isin(sub, te).any(), "training subset touches the test fold"
                f = (Y[sub]*S[sub]).sum(0)/np.maximum(S[sub].sum(0), 1)
                assert not np.isnan(f).any()
                pf[te] = f[None]
                pm[te] = fit_predict(T, M, Y, S, sub, te, dev)
            mm = per_complex(Y, pm, S); ff = per_complex(Y, pf, S)
            store[ds] = (pm, pf, mm, ff)
            d, lo, hi, p = cluster_bootstrap_diff(mm["auprc"], ff["auprc"], clu,
                                                  n_boot=N_BOOT, seed=BOOT_SEED)
            draw_rows.append(dict(threshold=thr, draw_seed=ds, n=int(np.sum(~np.isnan(mm["auprc"]))),
                                  model=round(float(np.nanmean(mm["auprc"])), 4),
                                  freq=round(float(np.nanmean(ff["auprc"])), 4),
                                  diff=round(float(d), 4), lo=round(float(lo), 4),
                                  hi=round(float(hi), 4), p=round(float(p), 4)))
            print(f"  thr {thr}  draw {ds}: model={draw_rows[-1]['model']:.4f} "
                  f"freq={draw_rows[-1]['freq']:.4f} diff={d:+.4f} [{lo:+.4f},{hi:+.4f}] "
                  f"p={p:.4f}  ({time.time()-t0:.1f}s)", flush=True)
            for i in range(n):
                per_rows.append(dict(threshold=thr, draw_seed=ds, pdb=pdbs[i], clu70=clu[i],
                                     model_auprc=mm["auprc"][i], freq_auprc=ff["auprc"][i],
                                     n_pos=mm["n_pos"][i]))
        # mean-over-draws predictions, bootstrapped once
        pm = np.nanmean([store[ds][0] for ds in DRAW_SEEDS], 0)
        pf = np.nanmean([store[ds][1] for ds in DRAW_SEEDS], 0)
        mm = per_complex(Y, pm, S); ff = per_complex(Y, pf, S)
        d, lo, hi, p = cluster_bootstrap_diff(mm["auprc"], ff["auprc"], clu,
                                              n_boot=N_BOOT, seed=BOOT_SEED)
        draw_rows.append(dict(threshold=thr, draw_seed="mean_over_draws",
                              n=int(np.sum(~np.isnan(mm["auprc"]))),
                              model=round(float(np.nanmean(mm["auprc"])), 4),
                              freq=round(float(np.nanmean(ff["auprc"])), 4),
                              diff=round(float(d), 4), lo=round(float(lo), 4),
                              hi=round(float(hi), 4), p=round(float(p), 4)))
        print(f"thr {thr}  MEAN OVER DRAWS: diff={d:+.4f} [{lo:+.4f},{hi:+.4f}] p={p:.4f}", flush=True)
        np.savez_compressed(OUT/f"predictions_thr{thr}.npz",
                            pred_model=pm, pred_freq=pf, Y=Y, S=S, pdb=pdbs)

    pd.DataFrame(draw_rows).to_csv(OUT/"per_draw_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(per_rows).to_csv(OUT/"per_complex_metrics.tsv", sep="\t", index=False)

    # the comparison the plan declares: Delta_sim - Delta_rand, drawn against the frozen sweep
    sweep = pd.read_csv(SWEEP_CSV).set_index("threshold")
    dd = pd.DataFrame(draw_rows)
    summ = []
    for thr in THRESHOLDS:
        g = dd[(dd.threshold == thr) & (dd.draw_seed != "mean_over_draws")]
        mo = dd[(dd.threshold == thr) & (dd.draw_seed == "mean_over_draws")].iloc[0]
        s = sweep.loc[thr]
        summ.append(dict(threshold=thr,
                         sweep_model=s.model, sweep_freq=s.freq, delta_sim=s["diff"],
                         sweep_lo=s.lo, sweep_hi=s.hi, sweep_p=s.p,
                         rand_model_mean=round(g.model.mean(), 4), rand_freq_mean=round(g.freq.mean(), 4),
                         delta_rand_mean=round(g["diff"].mean(), 4),
                         delta_rand_min=round(g["diff"].min(), 4), delta_rand_max=round(g["diff"].max(), 4),
                         delta_rand_range=round(g["diff"].max()-g["diff"].min(), 4),
                         delta_rand_meanpred=mo["diff"], meanpred_lo=mo.lo, meanpred_hi=mo.hi,
                         meanpred_p=mo.p,
                         delta_sim_minus_delta_rand=round(s["diff"] - g["diff"].mean(), 4)))
    sm = pd.DataFrame(summ)
    sm.to_csv(OUT/"size_matched_summary.tsv", sep="\t", index=False)
    print("\n=== Delta_sim (similarity-restricted) vs Delta_rand (size-matched random) ===")
    print(sm.to_string(index=False))
    unres = pd.read_csv(SWEEP_CSV).set_index("threshold").loc[1.01]
    print(f"\nunrestricted reference (threshold 1.01, 139-140 training complexes): "
          f"model={unres.model} freq={unres.freq} diff={unres['diff']}")
    print(f"\ntotal {time.time()-t_all:.1f}s")


if __name__ == "__main__":
    main()
