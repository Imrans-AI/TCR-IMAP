#!/usr/bin/env python3
"""The `interpretable` flag (criteria document, '"Very few" fixed at 10', 24 October 2026):
interpretable = rho defined and n_nonzero_occupancy >= 10. A description rule only - nothing is suppressed and no
other column changes. Synthetic inputs with known answers, directly and through statistics_phase, with the 9/10
boundary placed inside real strata. No trajectory and no model score is read.

usage: python md/analysis/tests/test_interpretable.py
"""
import sys
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import stats as st                       # noqa: E402
from pipeline import statistics_phase    # noqa: E402

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' - ' + detail) if detail else ''}")


print("=== add_interpretable on a hand-made table ===")
df = pd.DataFrame(dict(rho=[0.3, 0.3, 0.3, np.nan, np.nan, -0.2],
                       perm_p_approximate=[0.01, 0.01, 0.4, np.nan, np.nan, 0.5],
                       n_nonzero_occupancy=[9, 10, 400, 0, np.nan, 11]))
out = st.add_interpretable(df)
check("threshold is the predeclared 10", st.MIN_NONZERO_INTERPRETABLE == 10)
check("9 non-zero -> False, 10 -> True (boundary: fewer than 10 is uninterpretable)",
      out.interpretable.tolist()[:2] == [False, True])
check("many non-zero -> True; negative rho is no different", bool(out.interpretable[2]) and bool(out.interpretable[5]))
check("undefined rho -> False; no counts at all (NaN) -> False", out.interpretable.tolist()[3:5] == [False, False])
check("flag is boolean", out.interpretable.dtype == bool)
check("no other column changed, row count unchanged",
      out.drop(columns="interpretable").equals(df) and len(out) == len(df))

print("\n=== through statistics_phase ===")
rng = np.random.default_rng(9)
REPS = (1, 2, 3)
labels, occ, dists = {}, {}, {}
N_MID = {"A": 9, "B": 10}                    # non-zero occupancies placed in each system's 6-8 A stratum
for pdb in ("A", "B"):
    Y = np.zeros((40, 10), int); Y[:2, :2] = 1
    S = np.ones((40, 10), int)
    labels[pdb] = (Y, S)
    d = np.full((40, 10), 12.0); d[Y.astype(bool)] = 3.0
    d[2, :] = 5.0                             # 4.5-6 A: 10 pairs  -> untested (< 20)
    d[3:9, :] = 7.0                           # 6-8 A: 60 pairs
    dists[pdb] = d
    mid = np.argwhere(d == 7.0)
    pick = mid[rng.choice(len(mid), N_MID[pdb], replace=False)]
    for r in REPS:
        o = np.zeros((40, 10)); o[Y.astype(bool)] = 0.9
        o[2, :] = rng.random(10)
        for i, j in pick:
            o[i, j] = 0.05 + 0.1 * rng.random()          # non-zero in EVERY replicate -> counted in the mean
        # > 8 A stays all zero -> undefined
        n_dock = 5000 if pdb == "B" else 10                # A never contributes to the docked analysis
        occ[(pdb, r)] = dict(occ=o, occ_docked=o, n_docked=n_dock, occ_docked_identity=o, n_docked_identity=n_dock)
scores = {"TCR-IMAP": {p: rng.random((40, 10)) for p in labels},
          "pair-frequency prior": {p: rng.random((40, 10)) for p in labels}}
res = statistics_phase(occ, labels, scores, ("A", "B"), REPS, n_perm=199, dists=dists)
for name in ("rho", "strata", "docked_rho"):
    t = res[name]
    check(f"{name}: has `interpretable`, and it equals (rho defined) & (n_nonzero >= 10) on every row",
          "interpretable" in t.columns and
          (t.interpretable == (t.rho.notna() & (t.n_nonzero_occupancy.astype(float) >= 10))).all())
    check(f"{name}: no row suppressed, every rho and p still reported where defined",
          t.rho.notna().sum() >= (t.interpretable).sum())
s = res["strata"]
a_mid = s[(s.pdb == "A") & (s.stratum == "6-8 A")]
b_mid = s[(s.pdb == "B") & (s.stratum == "6-8 A")]
check("A 6-8 A: 9 non-zero -> rho and p reported, interpretable False",
      (a_mid.n_nonzero_occupancy == 9).all() and a_mid.rho.notna().all() and a_mid.perm_p_approximate.notna().all()
      and (~a_mid.interpretable).all())
check("B 6-8 A: 10 non-zero -> rho and p reported, interpretable True",
      (b_mid.n_nonzero_occupancy == 10).all() and b_mid.rho.notna().all() and b_mid.interpretable.all())
near = s[(s.level == "system") & (s.stratum == "4.5-6 A")]
check("4.5-6 A with 10 pairs: still 'not tested' (existing rule unchanged), interpretable False",
      near.note.str.contains("not tested").all() and near.rho.isna().all() and (~near.interpretable).all())
far = s[(s.level == "system") & (s.stratum == "> 8 A")]
check("> 8 A all zero: still undefined NaN/NaN (existing rule unchanged), interpretable False",
      far.rho.isna().all() and far.perm_p_approximate.isna().all() and (~far.interpretable).all())
pooled_mid = s[(s.level == "pooled") & (s.stratum == "6-8 A")]
check("pooled 6-8 A: 9 + 10 = 19 non-zero -> interpretable True",
      (pooled_mid.n_nonzero_occupancy == 19).all() and pooled_mid.interpretable.all())
dr = res["docked_rho"]
check("docked, A (no contributing replicate): NaN rho and counts, interpretable False",
      (~dr[dr.pdb == "A"].interpretable).all() and dr[dr.pdb == "A"].rho.isna().all())
check("docked, B: 10 non-zero mid pairs + 10 near pairs -> interpretable True",
      dr[dr.pdb == "B"].interpretable.all())
check("the flag is the LAST column of each table (no existing column moved)",
      all(res[n].columns[-1] == "interpretable" for n in ("rho", "strata", "docked_rho")))
no_d = statistics_phase(occ, labels, scores, ("A", "B"), REPS, n_perm=99)
check("called without distances (as the smoke test does): strata table empty, no error",
      no_d["strata"].empty and "interpretable" in no_d["rho"].columns)

print(f"\n{len(PASS)}/{len(PASS) + len(FAIL)} checks passed")
if FAIL:
    print("FAILED: " + ", ".join(FAIL))
sys.exit(1 if FAIL else 0)
