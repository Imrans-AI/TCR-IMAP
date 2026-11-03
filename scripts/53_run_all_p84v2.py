#!/usr/bin/env python3
"""p84v2 runner for scripts/24b_run_all_p82v1.py (pairwise TCR-IMAP vs pair-frequency prior, five protocols; PairNet from 19_train_pairmap). Executes that script's own source with the p84v2 data
layout and position-placed tokens; see scripts/53_runner_common.py for the exact, counted edits.
Default --check: builds the consumer's inputs and instantiates its models, no training.
--run trains; NOT used by the rebuild task (results/review/CORRECTION_PLAN_P84V2.md).
usage: python scripts/53_run_all_p84v2.py [--check | --run]
"""
import sys
sys.path.insert(0, "scripts")
from importlib import import_module
R = import_module("53_runner_common")

def check(ns):
    Y, S, T, M, meta = ns["load"]()
    return dict(T=T, M=M, pdbs=meta.pdb, arrays=dict(Y=Y, S=S))

sys.exit(0 if R.run("scripts/24b_run_all_p82v1.py", check_entry=check) in (True, None) else 1)
