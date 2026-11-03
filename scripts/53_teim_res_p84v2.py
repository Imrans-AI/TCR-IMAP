#!/usr/bin/env python3
"""p84v2 runner for scripts/36b_teim_res_p82v1.py (TEIM-Res comparator, pairwise, tcr70). Executes that script's own source with the p84v2 data
layout and position-placed tokens; see scripts/53_runner_common.py for the exact, counted edits.
Default --check: builds the consumer's inputs and instantiates its models, no training.
--run trains; NOT used by the rebuild task (results/review/CORRECTION_PLAN_P84V2.md).
usage: python scripts/53_teim_res_p84v2.py [--check | --run]

Determinism (added 23 Oct 2026, results/review/TEIM_DETERMINISM_REPORT.md). TEIMRes is the only comparator with 2-D
convolutions; under torch defaults cuDNN picks a nondeterministic weight-gradient kernel for its 3x3 dilated
Conv2d layers, so two runs with the same seeds gave different predictions (macro-AUPRC 0.2618 vs 0.2697). This
runner therefore sets, before 36b's code runs: CUBLAS_WORKSPACE_CONFIG=:4096:8 (before torch is imported),
torch.backends.cudnn.deterministic = True, torch.backends.cudnn.benchmark = False and
torch.use_deterministic_algorithms(True, warn_only=True). Any operation PyTorch reports as having no deterministic
implementation is printed once as "DETERMINISM-WARNING: ..." in the stage log. 36b itself is not edited.
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

def check(ns):
    L = ns["main"]()
    return dict(T=L["T"], M=L["M"], pdbs=L["meta"].pdb, arrays=dict(Y=L["Y"], S=L["S"]))

ok = R.run("scripts/36b_teim_res_p82v1.py", extra=[("mhc_len=Y.shape[2]; assert mhc_len == 82", "mhc_len=Y.shape[2]; assert mhc_len == 84", 1)],
           hook_after=r"^\s*M = _P84\.M_for\(meta\.pdb\)", check_entry=check)
print(f"determinism warnings raised: {len(_seen)}" + ("" if not _seen else " (see DETERMINISM-WARNING lines above)"), flush=True)
sys.exit(0 if ok in (True, None) else 1)
