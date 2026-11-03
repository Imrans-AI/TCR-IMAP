#!/usr/bin/env python3
"""n_nonzero_occupancy and n_occupancy_ge_025 beside every correlation (criteria document, "Occupancy counts
reported beside every correlation", 24 October 2026). Synthetic inputs with known answers: the helper itself
(including the 0.25 boundary and a tiny positive value), then every row level of the rho, strata and docked-rho
tables exactly as statistics_phase produces them, with counts recomputed independently from the inputs.
No trajectory and no model score is read.

usage: python md/analysis/tests/test_occupancy_counts.py
"""
import sys
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import stats as st                       # noqa: E402
from pipeline import statistics_phase    # noqa: E402

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' - ' + detail) if detail else ''}")


print("=== the helper ===")
c = st.occupancy_counts([0.0, 0.0, 1e-6, 0.2499, 0.25, 0.9, 0.9])
check("known group: distinct 5, non-zero 5, >= 0.25 3", c == dict(n_distinct_occupancy=5, n_nonzero_occupancy=5,
                                                                    n_occupancy_ge_025=3), str(c))
check("0.25 itself counts as >= 0.25 (boundary inclusive), 0.2499 does not",
      st.occupancy_counts([0.25])["n_occupancy_ge_025"] == 1 and st.occupancy_counts([0.2499])["n_occupancy_ge_025"] == 0)
check("a tiny positive occupancy (one frame in 9,001) counts as non-zero",
      st.occupancy_counts([1 / 9001])["n_nonzero_occupancy"] == 1)
check("all-zero group: 1 distinct, 0 non-zero, 0 >= 0.25",
      st.occupancy_counts(np.zeros(500)) == dict(n_distinct_occupancy=1, n_nonzero_occupancy=0, n_occupancy_ge_025=0))
check("empty group: all counts 0", st.occupancy_counts([]) == dict(n_distinct_occupancy=0, n_nonzero_occupancy=0,
                                                                    n_occupancy_ge_025=0))
check("the >= 0.25 threshold is the fixed MEANINGFUL constant", st.MEANINGFUL == 0.25)

print("\n=== through statistics_phase, every table and row level ===")
rng = np.random.default_rng(5)
REPS = (1, 2, 3)
labels, occ, dists = {}, {}, {}
for pdb in ("Z", "V"):
    Y = np.zeros((30, 10), int); Y[:2, :2] = 1
    S = np.ones((30, 10), int); S[29, :] = 0                    # one unresolved row, excluded everywhere
    labels[pdb] = (Y, S)
    d = np.full((30, 10), 12.0); d[Y.astype(bool)] = 3.0; d[2:6, :] = 5.0; d[6:12, :] = 7.0
    dists[pdb] = d
    for r in REPS:
        o = np.zeros((30, 10)); o[Y.astype(bool)] = 0.9
        o[2:6, :] = rng.random((4, 10)) * 0.6                    # near stratum: some >= 0.25
        o[6:12, :] = np.where(rng.random((6, 10)) < 0.2, 0.05, 0.0)   # 6-8 A: sparse, all < 0.25
        if pdb == "V":
            o[12:, :] = np.where(rng.random((18, 10)) < 0.05, 0.3, 0.0)   # far: a few >= 0.25
        o[29, :] = 0.99                                          # unresolved row: must never be counted
        nd = int(rng.integers(1000, 9001))
        occ[(pdb, r)] = dict(occ=o, occ_docked=o * 0.5, n_docked=nd if pdb == "V" else 10,   # Z never contributes
                             occ_docked_identity=o * 0.5, n_docked_identity=nd if pdb == "V" else 10)
scores = {"TCR-IMAP": {p: rng.random((30, 10)) for p in labels},
          "pair-frequency prior": {p: rng.random((30, 10)) for p in labels}}
res = statistics_phase(occ, labels, scores, ("Z", "V"), REPS, n_perm=199, dists=dists)
COLS = ["n_distinct_occupancy", "n_nonzero_occupancy", "n_occupancy_ge_025"]


def expect(vals):
    v = np.asarray(vals, float)
    return [len(np.unique(v)), int((v > 0).sum()), int((v >= 0.25).sum())]


neg = {p: st.neg_mask(*labels[p]) for p in labels}
sys_occ = {p: np.mean([occ[(p, r)]["occ"] for r in REPS], axis=0) for p in labels}
for name in ("rho", "strata", "docked_rho"):
    check(f"{name} table carries all three count columns", all(c in res[name].columns for c in COLS))

ok = True
for _, row in res["rho"].iterrows():
    if row.level == "system":
        want = expect(sys_occ[row.pdb][neg[row.pdb]])
    elif row.level == "replicate":
        want = expect(occ[(row.pdb, int(row.rep))]["occ"][neg[row.pdb]])
    else:
        want = expect(np.concatenate([sys_occ[p][neg[p]] for p in labels]))
    ok &= [int(row[c]) for c in COLS] == want
check("rho table: counts equal an independent recount at system, replicate and pooled level", ok)

ok = True
for _, row in res["strata"].iterrows():
    lo, hi = {n: (a, b) for n, a, b in st.DIST_STRATA}[row.stratum]
    ps = [row.pdb] if row.level == "system" else list(labels)
    vals = np.concatenate([sys_occ[p][neg[p] & (dists[p] >= lo) & (dists[p] < hi)] for p in ps])
    ok &= [int(row[c]) for c in COLS] == expect(vals)
check("strata table: counts equal an independent recount, per system and pooled", ok)
far_z = res["strata"][(res["strata"].pdb == "Z") & (res["strata"].stratum == "> 8 A")]
check("Z > 8 A: all-zero -> 0 non-zero, 0 >= 0.25, rho undefined",
      (far_z.n_nonzero_occupancy == 0).all() and (far_z.n_occupancy_ge_025 == 0).all() and far_z.rho.isna().all())
mid = res["strata"][(res["strata"].stratum == "6-8 A") & (res["strata"].level == "system")]
check("6-8 A: sparse non-zero occupancies all below 0.25 -> n_occupancy_ge_025 = 0 while n_nonzero > 0",
      (mid.n_occupancy_ge_025 == 0).all() and (mid.n_nonzero_occupancy > 0).all())

dr = res["docked_rho"]
z = dr[(dr.pdb == "Z")]
check("docked: Z has no contributing replicate -> all three counts NaN (no docked occupancy to count)",
      z[COLS].isna().all().all() and z.rho.isna().all())
v = dr[(dr.pdb == "V")]
want_v = expect(np.mean([occ[("V", r)]["occ_docked"] for r in REPS], axis=0)[neg["V"]])
check("docked: V counts equal a recount of its docked occupancy (both readings, both scores)",
      all([int(r[c]) for c in COLS] == want_v for _, r in v.iterrows()))
pooled = dr[dr.level == "pooled"]
check("docked pooled: counts over the contributing systems only (V)",
      all([int(r[c]) for c in COLS] == want_v for _, r in pooled.iterrows()))
check("no count gates reporting: every rho row is present whatever its counts",
      len(res["rho"]) == 2 * (2 * (1 + 3) + 1) and len(res["strata"]) == 2 * (2 * 3 + 3))
check("the unresolved row (occupancy 0.99) is never counted",
      res["rho"][res["rho"].level == "system"].n_crystal_negative.eq(int(neg["Z"].sum())).all())

print(f"\n{len(PASS)}/{len(PASS) + len(FAIL)} checks passed")
if FAIL:
    print("FAILED: " + ", ".join(FAIL))
sys.exit(1 if FAIL else 0)
