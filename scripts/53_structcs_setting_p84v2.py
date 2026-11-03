#!/usr/bin/env python3
"""p84v2 runner for scripts/31b_structcs_setting_p82v1.py (comparator-subset (88 structures) at 4.5 and 6 A: Ours, Struct-TCS, Seq-TCS). Executes that script's own source with the p84v2 data
layout and position-placed tokens; see scripts/53_runner_common.py for the exact, counted edits.
Default --check: builds the consumer's inputs and instantiates its models, no training.
--run trains; NOT used by the rebuild task (results/review/CORRECTION_PLAN_P84V2.md).
usage: python scripts/53_structcs_setting_p84v2.py [--check | --run]

Cluster-disjoint validation split (added 24 Oct 2026; predeclared in CORRECTION_PLAN_P84V2.md, section
"Comparator-subset validation split: the same cluster-disjoint correction"). 31b selects its three per-residue
models on a random 15 % of the non-test complexes, which is not disjoint by clu70 cluster - the same defect
corrected in the benchmark runner in Part B, and the reason Figure 4F would otherwise report a different
model-selection methodology from panels A, B and E of the same figure. This runner replaces that draw with whole
clu70 clusters, drawn at random from the fold's training complexes until at least 15 % are held out, and asserts
cluster-disjointness for every fold. The test folds themselves are drawn from the generator before the validation
loop and are therefore unchanged. The same code path serves all three methods.

It also persists the per-complex AUPRC arrays that 31b computes and discards, so the comparison can be made at the
paper's 20,000-resample seed-42 standard rather than only at the 2,000-resample values recorded at run time. That
records more and computes nothing differently.

BENCH_VALIDATION_SPLIT=random keeps 31b's original random draw, under this runner's determinism settings and p84v2
data, to reproduce the superseded values for the SI sensitivity analysis. It is NOT the reported configuration.

Determinism: as established for TEIM-Res (results/review/TEIM_DETERMINISM_REPORT.md).
"""
import os
os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"   # must precede the first cuBLAS call, so before torch is imported
import sys
import warnings
import torch
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
torch.use_deterministic_algorithms(True, warn_only=True)
_seen, _show = set(), warnings.showwarning
def _det_warning(message, category, filename, lineno, file=None, line=None):
    txt = str(message)
    if "determinis" in txt:
        if txt not in _seen:
            _seen.add(txt); print(f"DETERMINISM-WARNING: {category.__name__} at {filename}:{lineno}: {txt}", flush=True)
    else:
        _show(message, category, filename, lineno, file, line)
warnings.showwarning = _det_warning
print(f"determinism: CUBLAS_WORKSPACE_CONFIG={os.environ['CUBLAS_WORKSPACE_CONFIG']} cudnn.deterministic={torch.backends.cudnn.deterministic} "
      f"cudnn.benchmark={torch.backends.cudnn.benchmark} deterministic_algorithms={torch.are_deterministic_algorithms_enabled()} "
      f"(warn_only={torch.is_deterministic_algorithms_warn_only_enabled()}); torch {torch.__version__}, cuDNN {torch.backends.cudnn.version()}",
      flush=True)

sys.path.insert(0, "scripts")
from importlib import import_module
R = import_module("53_runner_common")

# Inserted before run_setting in the patched source; `np` is the consumer's own import. Identical in behaviour to
# the helper used by the benchmark runner (scripts/53_benchmark_p84v2.py).
VALIDATION_HELPER = '''def cluster_disjoint_validation(clu, rest, rng, frac=0.15):
    """Validation complexes as whole clu70 clusters (p84v2 correction, predeclared in
    results/review/CORRECTION_PLAN_P84V2.md). Clusters present among the training complexes `rest` are drawn at
    random without replacement until at least `frac` of those complexes are held out; every complex of a drawn
    cluster is held out, so no validation complex shares a cluster with any remaining training complex.
    One draw per fold from the caller's generator."""
    cl = clu[rest]
    uc = np.unique(cl)
    target = frac * len(rest)
    sel, n = [], 0
    for c in uc[rng.permutation(len(uc))]:
        if n >= target:
            break
        sel.append(c)
        n += int((cl == c).sum())
    return rest[np.isin(cl, sel)]


'''

OLD_SPLIT = """        vm=rng.random(len(rest))<0.15; tr,va=rest[~vm],rest[vm]
        if len(va)<3 or len(tr)<15: continue
"""
NEW_SPLIT = """        va = cluster_disjoint_validation(meta.clu70.values, rest, rng, 0.15)
        tr = np.setdiff1d(rest, va)
        if len(va)<3 or len(tr)<15: continue
        _tc, _vc = set(meta.clu70.values[tr]), set(meta.clu70.values[va])
        assert not (_tc & _vc), f"validation not cluster-disjoint ({setting}): shared clu70 {sorted(_tc & _vc)}"
        print(f"   [{setting}] fold: test {len(te)}, train {len(tr)}, validation {len(va)} complexes in "
              f"{len(_vc)} clu70 clusters = {100*len(va)/len(rest):.1f}% of {len(rest)} training complexes", flush=True)
"""
SENSITIVITY_SPLIT = """        vm=rng.random(len(rest))<0.15; tr,va=rest[~vm],rest[vm]
        if len(va)<3 or len(tr)<15: continue
        print(f"   [{setting}] SENSITIVITY random split: test {len(te)}, train {len(tr)}, validation {len(va)} "
              f"complexes = {100*len(va)/len(rest):.1f}% of {len(rest)} training complexes, NOT cluster-disjoint "
              f"({len(set(meta.clu70.values[tr]) & set(meta.clu70.values[va]))} clu70 clusters on both sides)", flush=True)
"""

# Per-seed predictions (added 24 Oct 2026 for the Figure 4F stability check). The three seeds are fitted in
# exactly the same order with exactly the same arguments as before and the SAME list is averaged, so the
# averaged output is unchanged; the per-seed matrices are merely kept instead of being discarded.
OLD_SEEDS = """        for name,(cls,kw) in MODELS.items():
            P[name][te]=np.mean([fit(cls,kw,T,M,A,Y,S,tr,va,te,s,dev) for s in SEEDS],0)
"""
NEW_SEEDS = """        for name,(cls,kw) in MODELS.items():
            _preds=[fit(cls,kw,T,M,A,Y,S,tr,va,te,s,dev) for s in SEEDS]
            for _s,_pr in zip(SEEDS,_preds): _PSEED.setdefault(name,{}).setdefault(_s,np.full(Y.shape,np.nan,np.float32))[te]=_pr
            P[name][te]=np.mean(_preds,0)
"""
OLD_PSINIT = """    P={k:np.full(Y.shape,np.nan,np.float32) for k in MODELS}
    used=np.zeros(n,bool)
"""
NEW_PSINIT = """    P={k:np.full(Y.shape,np.nan,np.float32) for k in MODELS}
    _PSEED={}
    used=np.zeros(n,bool)
"""

# persist the per-complex arrays 31b computes and discards (no change to any computed value)
OLD_PC = """    mc={k:per_complex(Y[i],P[k][i],S[i]) for k in MODELS}
    line = dict(setting=setting, n_test=len(i))
"""
NEW_PC = """    mc={k:per_complex(Y[i],P[k][i],S[i]) for k in MODELS}
    globals()["_PC"] = {"pdb": [str(x) for x in meta.pdb.values[i]],
                        "clu70": [int(x) for x in meta.clu70.values[i]],
                        **{f"{k}_auprc": [float(v) for v in mc[k]["auprc"]] for k in MODELS}}
    _ms={(k,s):per_complex(Y[i],_PSEED[k][s][i],S[i]) for k in MODELS for s in SEEDS}
    globals()["_PSEEDPC"] = {f"{k}_seed{s}_auprc": [float(v) for v in _ms[(k,s)]["auprc"]]
                             for k in MODELS for s in SEEDS}
    import os as _os
    _d = "results/tables/_partial/p84v2/31b_perseed"; _os.makedirs(_d, exist_ok=True)
    np.savez(f"{_d}/{_UNIT}", pdb=meta.pdb.values[i], clu70=meta.clu70.values[i],
             **{f"{k}_seed{s}": _PSEED[k][s][i] for k in MODELS for s in SEEDS})
    line = dict(setting=setting, n_test=len(i))
"""
OLD_CACHE = """    for setting in ("their_random","our_tcr70"):
        t0=time.time()
        def compute():
            line = run_setting(setting, T,M,A,Y,S,meta,dev)
            return {"__json__": np.array(json.dumps(line))}
"""
NEW_CACHE = """    for setting in ("their_random","our_tcr70"):
        t0=time.time()
        globals()["_UNIT"] = f"{tag}_{setting}"
        def compute():
            line = run_setting(setting, T,M,A,Y,S,meta,dev)
            return {"__json__": np.array(json.dumps(line)),
                    "__percomplex__": np.array(json.dumps(_PC)),
                    "__perseed__": np.array(json.dumps(_PSEEDPC))}
"""

SPLIT = os.environ.get("BENCH_VALIDATION_SPLIT", "cluster")
assert SPLIT in ("cluster", "random"), f"BENCH_VALIDATION_SPLIT must be 'cluster' or 'random', got {SPLIT!r}"
print(f"validation split: {SPLIT}" + ("" if SPLIT == "cluster" else
      "   <-- SENSITIVITY configuration, not the reported one"), flush=True)


def check(ns):
    T, M, A, Y, S, meta, n = ns["load"]("data/processed/pair_maps_p84v2.npz")
    T6, M6, A6, Y6, S6, meta6, n6 = ns["load"]("data/processed/pair_maps_p84v2_6A.npz")
    assert (T == T6).all() and (M == M6).all() and list(meta.pdb) == list(meta6.pdb)
    return dict(T=T, M=M, pdbs=meta.pdb, arrays=dict(Y_per_residue_4p5A=Y, Y_per_residue_6A=Y6, A_graph=A))


ok = R.run("scripts/31b_structcs_setting_p82v1.py",
           extra=[("def run_setting(", VALIDATION_HELPER + "def run_setting(", 1),
                  (OLD_SPLIT, NEW_SPLIT if SPLIT == "cluster" else SENSITIVITY_SPLIT, 1),
                  (OLD_PSINIT, NEW_PSINIT, 1),
                  (OLD_SEEDS, NEW_SEEDS, 1),
                  (OLD_PC, NEW_PC, 1),
                  (OLD_CACHE, NEW_CACHE, 1)],
           check_entry=check)
print(f"determinism warnings raised: {len(_seen)}" + ("" if not _seen else " (see DETERMINISM-WARNING lines above)"), flush=True)
sys.exit(0 if ok in (True, None) else 1)
