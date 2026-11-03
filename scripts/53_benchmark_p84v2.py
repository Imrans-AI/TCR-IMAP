#!/usr/bin/env python3
"""p84v2 runner for scripts/28b_benchmark_p82v1.py (per-residue benchmark: Ours, Struct-TCS, Seq-TCS (21_struct_tcs), position prior, NN transfer). Executes that script's own source with the p84v2 data
layout and position-placed tokens; see scripts/53_runner_common.py for the exact, counted edits.
Default --check: builds the consumer's inputs and instantiates its models, no training.
--run trains; NOT used by the rebuild task (results/review/CORRECTION_PLAN_P84V2.md).
The NN-transfer baseline uses no tokens: it copies the label rows of the training complex with the
highest Levenshtein ratio on the frozen complex_metadata.tcr_seq, so under p84v2 it copies p84v2 rows.
usage: python scripts/53_benchmark_p84v2.py [--check | --run]

Cluster-disjoint validation split (added 24 Oct 2026; predeclared in CORRECTION_PLAN_P84V2.md, section
"Per-residue validation split: cluster-disjoint correction"). 28b early-stops the three per-residue models on a
random 15 % of the non-test complexes, which is not disjoint by clu70 cluster, so a validation complex can be a
near-duplicate of a complex the model trains on. This runner replaces that draw with whole clu70 clusters, drawn at
random from the fold's training complexes until at least 15 % of them are held out, and asserts cluster-disjointness
for every fold. The same code path serves all three per-residue models, so the correction applies identically to
TCR-IMAP's per-residue variant, Struct-TCS and Seq-TCS. Nothing else changes: same RNG seed discipline
(one default_rng(0) per protocol, one draw per fold in fold order), same guards, epochs, patience, optimiser and
seeds. 28b itself is not edited.

Determinism: as established for TEIM-Res (results/review/TEIM_DETERMINISM_REPORT.md), any stage re-run under this
correction sets CUBLAS_WORKSPACE_CONFIG=:4096:8 (before torch is imported), cudnn.deterministic, cudnn.benchmark off
and torch.use_deterministic_algorithms(True, warn_only=True). Operations with no deterministic implementation (the
per-residue models use cuDNN LSTMs, unlike TEIM-Res) are printed once as "DETERMINISM-WARNING: ..." and counted at
exit, so the stage log records exactly what could not be made deterministic.
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

# Inserted before run_protocol in the patched source; `np` is the consumer's own import.
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

OLD_SPLIT = """        vmask = rng.random(len(rest))<0.15
        tr, va = rest[~vmask], rest[vmask]
        if len(va)<3 or len(tr)<20: continue
"""
# BENCH_VALIDATION_SPLIT=random keeps 28b's original random draw, under this runner's determinism settings and
# p84v2 data. It exists to reproduce the superseded random-split results for the SI sensitivity analysis the
# predeclaration requires, and to separate the effect of the split from the effect of the determinism settings.
# It is NOT the reported configuration; the default is the cluster-disjoint split.
SENSITIVITY_SPLIT = """        vmask = rng.random(len(rest))<0.15
        tr, va = rest[~vmask], rest[vmask]
        if len(va)<3 or len(tr)<20: continue
        print(f"   [{pname}] SENSITIVITY random split: test {len(te)}, train {len(tr)}, validation {len(va)} "
              f"complexes = {100*len(va)/len(rest):.1f}% of {len(rest)} training complexes, NOT cluster-disjoint "
              f"({len(set(meta.clu70.values[tr]) & set(meta.clu70.values[va]))} clu70 clusters on both sides)", flush=True)
"""
NEW_SPLIT = """        va = cluster_disjoint_validation(meta.clu70.values, rest, rng, 0.15)
        tr = np.setdiff1d(rest, va)
        if len(va)<3 or len(tr)<20: continue
        _tc, _vc = set(meta.clu70.values[tr]), set(meta.clu70.values[va])
        assert not (_tc & _vc), f"validation not cluster-disjoint ({pname}): shared clu70 {sorted(_tc & _vc)}"
        print(f"   [{pname}] fold: test {len(te)}, train {len(tr)}, validation {len(va)} complexes in "
              f"{len(_vc)} clu70 clusters = {100*len(va)/len(rest):.1f}% of {len(rest)} training complexes", flush=True)
"""

def check(ns):
    Y, S, T, M, A, meta, mhc_len = ns["load"]()
    return dict(T=T, M=M, pdbs=meta.pdb, arrays=dict(Y_per_residue=Y, S_per_residue=S, A_graph=A))

SPLIT = os.environ.get("BENCH_VALIDATION_SPLIT", "cluster")
assert SPLIT in ("cluster", "random"), f"BENCH_VALIDATION_SPLIT must be 'cluster' or 'random', got {SPLIT!r}"
print(f"validation split: {SPLIT}" + ("" if SPLIT == "cluster" else
      "   <-- SENSITIVITY configuration, not the reported one"), flush=True)

ok = R.run("scripts/28b_benchmark_p82v1.py",
           extra=[("assert mhc_len == 82", "assert mhc_len == 84", 1),
                  ("def run_protocol(", VALIDATION_HELPER + "def run_protocol(", 1),
                  (OLD_SPLIT, NEW_SPLIT if SPLIT == "cluster" else SENSITIVITY_SPLIT, 1)],
           check_entry=check)
print(f"determinism warnings raised: {len(_seen)}" + ("" if not _seen else " (see DETERMINISM-WARNING lines above)"), flush=True)
sys.exit(0 if ok in (True, None) else 1)
