#!/usr/bin/env python3
"""p84v2 runner for scripts/16b_docking_angle_p82v1.py (AngleNet docking-angle regressor (Fig 3C/D; a separate model, not TCR-IMAP)). Executes that script's own source with the p84v2 data
layout and position-placed tokens; see scripts/53_runner_common.py for the exact, counted edits.
Default --check: builds the consumer's inputs and instantiates its models, no training.
--run trains; NOT used by the rebuild task (results/review/CORRECTION_PLAN_P84V2.md).
AngleNet's CDR3 definition is kept EXPLICITLY: 16b read contact_sites_p82v1.parquet, whose CDR3
strings stop at IMGT 117 (11b), not 124. Tokens are placed by position on the p84v2 rows and the
CDR3 rows above IMGT 117 (FR4 118-124) are padding: ANGLENET_CDR3_LAST_IMGT = 117 below. Docking
angles are still read from contact_sites_p82v1.parquet (input only, unchanged).
usage: python scripts/53_docking_angle_p84v2.py [--check | --run]
"""
import sys
sys.path.insert(0, "scripts")
from importlib import import_module
R = import_module("53_runner_common")

ANGLENET_CDR3_LAST_IMGT = 117


def check(ns):
    L = ns["main"]()
    return dict(T=L["T"], M=L["M"], pdbs=L["df"].pdb, arrays=dict(y_angle=L["y"]) if "y" in L else {})

sys.exit(0 if R.run("scripts/16b_docking_angle_p82v1.py", cdr3_last_imgt=ANGLENET_CDR3_LAST_IMGT,
                    extra=[("    mhc_len = df.mhc_contact.str.len().iloc[0]\n", "    mhc_len = 84   # p84v2: 84 MHC columns (82 + 1061A + 1072A)\n", 1),
                           ("    assert (df.mhc_contact.str.len()==mhc_len).all() and mhc_len==82\n", "    assert mhc_len == 84\n", 1)],
                    hook_after=r"^\s*M = _P84\.M_for\(df\.pdb\)", check_entry=check) in (True, None) else 1)
