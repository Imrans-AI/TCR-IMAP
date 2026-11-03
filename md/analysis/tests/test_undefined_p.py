#!/usr/bin/env python3
"""Undefined correlations report undefined p (criteria document, "Undefined correlations are reported as
undefined", 24 October 2026). Synthetic inputs with known answers for every correlation the run computes:
constant occupancy, all-zero occupancy, n = 1, n = 2 and a normal case, for
  per-replicate rho (stats.rho), per-system rho + p (stats.perm_p), pooled rho + p (stats.pooled_perm_p),
  stratified rho + p (stats.rho_by_stratum, stats.pooled_rho_by_stratum),
and the same through pipeline.statistics_phase exactly as the run calls it, including the prior control and
both docked-frame readings. No trajectory and no model score is read.

usage: python md/analysis/tests/test_undefined_p.py
"""
import sys
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import stats as st                       # noqa: E402
from pipeline import statistics_phase    # noqa: E402

PASS, FAIL = [], []
NP = 999


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' - ' + detail) if detail else ''}")


def undefined(r, p):
    return np.isnan(r) and np.isnan(p)


rng = np.random.default_rng(3)


def group(n, kind):
    """n crystal-negative pairs in one column; returns score, occ, Y, S."""
    Y = np.zeros((n, 1), int); S = np.ones((n, 1), int)
    sc = rng.random((n, 1))
    oc = {"constant": np.full((n, 1), 0.37), "all-zero": np.zeros((n, 1)),
          "normal": rng.random((n, 1))}[kind]
    return sc, oc, Y, S


CASES = [("constant occupancy", 50, "constant"), ("all-zero occupancy", 50, "all-zero"),
         ("n = 1", 1, "normal"), ("n = 0", 0, "normal")]

print("=== undefined cases: rho and p both NaN, for every statistic ===")
for label, n, kind in CASES:
    sc, oc, Y, S = group(n, kind)
    check(f"{label}: per-replicate rho is NaN", np.isnan(st.rho(sc, oc, Y, S)))
    check(f"{label}: per-system rho and p are NaN", undefined(*st.perm_p(sc, oc, Y, S, n_perm=NP)))
    check(f"{label}: pooled rho and p are NaN", undefined(*st.pooled_perm_p([(sc, oc, Y, S)], n_perm=NP)))
    keep, st.MIN_STRATUM_N = st.MIN_STRATUM_N, 0          # reach the computation, not the n < 20 guard
    d = np.full((n, 1), 10.0)
    r = [x for x in st.rho_by_stratum(sc, oc, Y, S, d, n_perm=NP) if x["stratum"] == "> 8 A"][0]
    q = [x for x in st.pooled_rho_by_stratum([(sc, oc, Y, S, d)], n_perm=NP) if x["stratum"] == "> 8 A"][0]
    st.MIN_STRATUM_N = keep
    check(f"{label}: stratum rho and p are NaN, note gives the reason",
          undefined(r["rho"], r["perm_p_approximate"]) and r["note"].startswith("rho undefined"), repr(r["note"]))
    # an EMPTY pooled group (n = 0) is caught first by the existing "not tested" guard (n < MIN_STRATUM_N or no
    # pairs at all); either way rho and p must both be NaN and the note must say why
    check(f"{label}: pooled stratum rho and p are NaN",
          undefined(q["rho"], q["perm_p_approximate"]) and (q["note"].startswith("rho undefined")
                                                            or (n == 0 and q["note"].endswith("not tested"))),
          repr(q["note"]))
    check(f"{label}: n and n_distinct_occupancy reported", r["n"] == n and r["n_distinct_occupancy"] == min(n, 1)
          if kind == "normal" else r["n"] == n and r["n_distinct_occupancy"] == 1,
          f"n {r['n']}, distinct {r['n_distinct_occupancy']}")

sc, oc, Y, S = group(50, "normal"); sc[:] = 0.5
check("constant SCORES: per-system rho and p are NaN", undefined(*st.perm_p(sc, oc, Y, S, n_perm=NP)))
check("constant SCORES: reason is 'scores constant'", st.undefined_reason(sc.ravel(), oc.ravel()) == "scores constant")
check("_p never returns the floor for an undefined observed rho", np.isnan(st._p(np.nan, np.zeros(NP), "two-sided"))
      and np.isnan(st._p(np.nan, np.zeros(NP), "greater")))

print("\n=== defined cases: unchanged, and honest ===")
# n = 2 with distinct values: rho is exactly +1 or -1; every permutation is as extreme, so two-sided p = 1
for sgn in (+1, -1):
    Y = np.zeros((2, 1), int); S = np.ones((2, 1), int)
    sc = np.array([[0.2], [0.8]]); oc = np.array([[0.1], [0.9]]) if sgn > 0 else np.array([[0.9], [0.1]])
    r, p = st.perm_p(sc, oc, Y, S, n_perm=NP)
    check(f"n = 2, rho {sgn:+d}: rho = {sgn:+d} exactly and two-sided p = 1", r == sgn and p == 1.0, f"rho {r}, p {p}")
# normal: rho matches scipy, p a proper probability
sc, oc, Y, S = group(200, "normal")
r, p = st.perm_p(sc, oc, Y, S, n_perm=NP)
check("normal: rho equals scipy's Spearman", abs(r - spearmanr(sc.ravel(), oc.ravel()).correlation) < 1e-12)
check("normal: p is finite and in (0, 1]", np.isfinite(p) and 0 < p <= 1, f"p {p:.4f}")
check("normal: per-replicate rho equals perm_p's rho", abs(st.rho(sc, oc, Y, S) - r) < 1e-12)
# a defined, perfectly monotone association still reaches the floor - the floor is right when rho is defined
oc2 = sc.copy()
r2, p2 = st.perm_p(sc, oc2, Y, S, n_perm=NP)
check("defined perfect association: rho = 1 and p at the floor 1/(1+n_perm)", abs(r2 - 1) < 1e-12 and p2 == 1 / (NP + 1))
# pooling a constant block with a varying block is DEFINED (the union varies) and equals scipy on the union
a = group(40, "all-zero"); b = group(40, "normal")
rp, pp = st.pooled_perm_p([a, b], n_perm=NP)
ref = spearmanr(np.concatenate([a[0].ravel(), b[0].ravel()]), np.concatenate([a[1].ravel(), b[1].ravel()])).correlation
check("pooled: an all-zero system pooled with a varying one is defined and equals scipy on the union",
      abs(rp - ref) < 1e-12 and np.isfinite(pp))
check("pooled: two all-zero systems are undefined", undefined(*st.pooled_perm_p([a, group(40, "all-zero")], n_perm=NP)))

print("\n=== through pipeline.statistics_phase, as the run calls it ===")
# two systems on a 30 x 10 grid; system 'Z' has zero occupancy everywhere except its crystal contacts,
# system 'V' varies. Distances put 'Z''s far pairs in the > 8 A stratum.
REPS = (1, 2, 3)
labels, occ, dists = {}, {}, {}
for pdb in ("Z", "V"):
    Y = np.zeros((30, 10), int); Y[:2, :2] = 1
    S = np.ones((30, 10), int)
    labels[pdb] = (Y, S)
    d = np.full((30, 10), 12.0); d[Y.astype(bool)] = 3.0; d[2:6, :] = 5.0; d[6:12, :] = 7.0
    dists[pdb] = d
    for r in REPS:
        if pdb == "Z":
            o = np.zeros((30, 10)); o[Y.astype(bool)] = 0.9
            o[2:6, :] = rng.random((4, 10)) * 0.3; o[6:12, :] = rng.random((6, 10)) * 0.1   # near strata vary
        else:
            o = rng.random((30, 10))
        occ[(pdb, r)] = dict(occ=o, occ_docked=o, n_docked=5000, occ_docked_identity=o, n_docked_identity=5000)
scores = {"TCR-IMAP": {p: rng.random((30, 10)) for p in labels},
          "pair-frequency prior": {p: rng.random((30, 10)) for p in labels}}
res = statistics_phase(occ, labels, scores, ("Z", "V"), REPS, n_perm=NP, dists=dists)
rho, strata, drho = res["rho"], res["strata"], res["docked_rho"]
for t, name in ((rho, "rho"), (strata, "strata"), (drho, "docked rho")):
    check(f"{name} table: no p-value for an undefined rho", (t[t.rho.isna()].perm_p_approximate.isna()).all())
    check(f"{name} table: carries n_distinct_occupancy", "n_distinct_occupancy" in t.columns)
check("docked rho table: carries n_crystal_negative", "n_crystal_negative" in drho.columns)
far = strata[(strata.stratum == "> 8 A") & (strata.pdb == "Z")]
check("system Z, > 8 A, both scores: rho and p NaN, n_distinct_occupancy 1, reason in note",
      len(far) == 2 and far.rho.isna().all() and far.perm_p_approximate.isna().all()
      and (far.n_distinct_occupancy == 1).all() and far.note.str.contains("occupancy constant").all())
near = strata[(strata.stratum != "> 8 A") & (strata.level == "system")]
check("every other system stratum is defined, with a finite p", near.rho.notna().all() and near.perm_p_approximate.notna().all())
pooled_far = strata[(strata.stratum == "> 8 A") & (strata.level == "pooled")]
check("pooled > 8 A: defined, because system V varies there", pooled_far.rho.notna().all() and pooled_far.perm_p_approximate.notna().all())
check("whole-system, replicate and pooled rho all defined (Z varies across its strata)",
      rho.rho.notna().all() and rho[rho.level != "replicate"].perm_p_approximate.notna().all())
check("docked rows: both readings and both scores defined", drho.rho.notna().all() and set(drho.reading) == {"count", "identity"})
check("prior control handled identically (same rows, same NaN pattern)",
      (strata[strata.score == "TCR-IMAP"].rho.isna().values == strata[strata.score == "pair-frequency prior"].rho.isna().values).all())

print(f"\n{len(PASS)}/{len(PASS) + len(FAIL)} checks passed")
if FAIL:
    print("FAILED: " + ", ".join(FAIL))
sys.exit(1 if FAIL else 0)
