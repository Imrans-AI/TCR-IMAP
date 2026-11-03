#!/usr/bin/env python3
"""Step 4: every statistic of the MD occupancy analysis on synthetic inputs with known answers.
No real score is combined with real MD data here (no real data is read at all).
usage: python md/analysis/tests/test_stats.py      (or: pytest md/analysis/tests)
"""
import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import stats as st  # noqa: E402
from occupancy import is_docked, in_window  # noqa: E402
from contacts import PairMinDist  # noqa: E402

RNG = np.random.default_rng(0)


def grid(score_vals, occ_vals):
    """Embed 1-D synthetic pairs into 84 x 82 matrices: the first n cells are crystal-negative,
    every other cell is unresolved (S = 0) and filled with adversarial values."""
    n = len(score_vals)
    S = np.zeros((84, 82), int); Y = np.zeros((84, 82), int)
    sc = np.full((84, 82), 1e6); oc = np.full((84, 82), -5.0)
    idx = np.unravel_index(np.arange(n), (84, 82))
    S[idx] = 1; sc[idx] = score_vals; oc[idx] = occ_vals
    return sc, oc, Y, S


def test_rho_monotone_and_reversed():
    x = RNG.random(200)
    assert abs(st.rho(*grid(x, x ** 3)) - 1.0) < 1e-12
    assert abs(st.rho(*grid(x, -x)) + 1.0) < 1e-12


def test_rho_hand_computed_with_ties():
    # ranks: score 1..5, occ (1.5, 1.5, 3.5, 3.5, 5) -> rho = 9 / sqrt(10 * 9) = 3 / sqrt(10)
    r = st.rho(*grid([1, 2, 3, 4, 5], [0, 0, 0.5, 0.5, 1]))
    assert abs(r - 3 / np.sqrt(10)) < 1e-12


def test_only_crystal_negative_pairs_enter():
    sc, oc, Y, S = grid(np.arange(50.0), np.arange(50.0))
    base = st.rho(sc, oc, Y, S)
    # add crystal-positive pairs with perfectly anti-correlated values: must not change rho
    idx = np.unravel_index(np.arange(50, 70), (84, 82))
    S[idx] = 1; Y[idx] = 1; sc[idx] = np.arange(20.0); oc[idx] = -np.arange(20.0)
    assert abs(base - 1.0) < 1e-12 and st.rho(sc, oc, Y, S) == base


def test_system_occupancy_is_replicate_mean():
    reps = [np.full((84, 82), v) for v in (0.2, 0.4, 0.9)]
    assert np.allclose(st.system_occupancy(reps), 0.5)


def test_pooled_rho_is_rho_over_union():
    a = grid([1, 2, 3], [0.1, 0.2, 0.3]); b = grid([10, 20, 30], [0.0, 0.05, 0.9])
    from scipy.stats import spearmanr
    ref = spearmanr([1, 2, 3, 10, 20, 30], [0.1, 0.2, 0.3, 0.0, 0.05, 0.9]).correlation
    assert abs(st.pooled_rho([a, b]) - ref) < 1e-12


def test_perm_p_floor_for_perfect_association():
    x = RNG.random(30)
    r, p = st.perm_p(*grid(x, x), n_perm=10_000)
    assert r == 1.0 and p == 1 / 10_001


def test_perm_p_matches_exact_enumeration():
    s = np.array([0.1, 0.4, 0.35, 0.8, 0.6, 0.2]); o = np.array([0.0, 0.3, 0.1, 0.5, 0.2, 0.25])
    rs, ro = st._rank(s), st._rank(o)
    obs = st._rho_from_ranks(rs, ro)
    null = np.array([st._rho_from_ranks(np.array(q), ro) for q in itertools.permutations(rs)])
    exact = (np.abs(null) >= abs(obs) - 1e-12).mean()
    r, p = st.perm_p(*grid(s, o), n_perm=20_000, seed=1)
    assert abs(r - obs) < 1e-12 and abs(p - exact) < 0.01, (p, exact)
    r, pg = st.perm_p(*grid(s, o), n_perm=20_000, seed=1, alternative="greater")
    exact_g = (null >= obs - 1e-12).mean()
    assert abs(pg - exact_g) < 0.01, (pg, exact_g)


def test_perm_p_calibrated_under_null():
    ps = []
    for i in range(200):
        x, y = RNG.random(80), RNG.random(80)
        ps.append(st.perm_p(*grid(x, y), n_perm=999, seed=i)[1])
    ps = np.array(ps)
    assert 0.02 <= (ps < 0.05).mean() <= 0.09, (ps < 0.05).mean()
    assert 0.40 <= np.median(ps) <= 0.60, np.median(ps)


def test_pooled_permutation_is_within_system():
    # Within each system score and occupancy are independent, but system B has higher scores AND
    # higher occupancy, so the pooled rho is large (~0.75). The within-system permutation keeps that
    # between-system part in its null, so over 200 such datasets its p is calibrated (~5 % below
    # 0.05); a global permutation of the same data would call every one significant.
    ps, rs = [], []
    for i in range(200):
        a = grid(RNG.random(150), RNG.random(150) * 0.1)
        b = grid(RNG.random(150) + 1, RNG.random(150) * 0.1 + 0.5)
        r, p = st.pooled_perm_p([a, b], n_perm=499, seed=i)
        ps.append(p); rs.append(r)
    ps = np.array(ps)
    assert min(rs) > 0.6, min(rs)
    assert 0.02 <= (ps < 0.05).mean() <= 0.09, (ps < 0.05).mean()
    # global permutation on the last dataset, for contrast: p at its floor
    s = np.concatenate([a[0][st.neg_mask(a[2], a[3])], b[0][st.neg_mask(b[2], b[3])]])
    o = np.concatenate([a[1][st.neg_mask(a[2], a[3])], b[1][st.neg_mask(b[2], b[3])]])
    assert st.perm_p(*grid(s, o), n_perm=499)[1] == 1 / 500
    # and genuine within-system association is detected
    x1, x2 = RNG.random(300), RNG.random(300)
    r2, p2 = st.pooled_perm_p([grid(x1, x1 + RNG.normal(0, .2, 300)), grid(x2, x2 + RNG.normal(0, .2, 300))], n_perm=2000)
    assert p2 == 1 / 2001 and r2 > 0.7


def test_top_k_known_answer():
    S = np.zeros((84, 82), int); Y = np.zeros((84, 82), int)
    sc = np.zeros((84, 82)); oc = np.zeros((84, 82))
    S[:10, :10] = 1                       # 100 resolved pairs
    Y[0, :5] = 1                          # k = 5 crystal-positive (given the highest scores: excluded)
    sc[0, :5] = 100.0
    neg_top = [(1, 0), (2, 0), (3, 0), (4, 0), (5, 0)]
    for n, (i, j) in enumerate(neg_top):
        sc[i, j] = 50.0 - n
    for i, j in neg_top[:3]:
        oc[i, j] = 0.30                   # 3 of the top 5 >= 0.25 (one exactly 0.25 below)
    oc[neg_top[3]] = 0.2499
    oc[9, 9] = 0.25                       # one more meaningful negative outside the top k
    r = st.top_k_fraction(sc, oc, Y, S)
    assert r["k"] == 5 and r["frac_top_k"] == 0.6 and r["n_top_k_ge"] == 3
    assert abs(r["base_rate"] - 4 / 95) < 1e-12 and r["tie_at_k"] is False
    sc[6, 0] = sc[5, 0]                   # tie between the 5th and 6th negative
    assert st.top_k_fraction(sc, oc, Y, S)["tie_at_k"] is True


def test_prior_control_uses_identical_code():
    x = RNG.random(100); occ = x.copy()
    model = grid(x, occ); prior = grid(-x, occ)
    assert abs(st.rho(*model) - 1.0) < 1e-12 and abs(st.rho(*prior) + 1.0) < 1e-12


def test_docked_rules():
    assert is_docked(50, 100) and not is_docked(49, 100)          # >= 50 %, boundary included
    assert st.docked_contributes(1000) and not st.docked_contributes(999)
    occs = [np.full((2, 2), v) for v in (0.2, 0.6, 1.0)]
    m, n = st.docked_system_occupancy(occs, [1000, 999, 5000])
    assert n == 2 and np.allclose(m, 0.6)                          # mean of reps 1 and 3 only
    assert st.docked_system_occupancy(occs, [10, 20, 999]) == (None, 0)


def test_window_selects_9001_frames():
    t = np.arange(0, 200_001, 20.0)
    assert len(t) == 10_001 and sum(in_window(x, 20_000, 200_000) for x in t) == 9_001


def test_positive_summary():
    S = np.zeros((84, 82), int); Y = np.zeros((84, 82), int); oc = np.zeros((84, 82))
    S[0, :4] = Y[0, :4] = 1; oc[0, :4] = [1.0, 0.9, 0.2, 0.0]
    r = st.positive_occupancy_summary(oc, Y, S)
    assert r["n"] == 4 and r["median"] == 0.55 and r["n_ge_meaningful"] == 2 and r["min"] == 0.0


def test_contact_function_synthetic_trajectory_with_pbc():
    # slot 0: two heavy atoms; MHC position 0: three atoms; one atom pair approaches to 4.0 A in
    # 3 of 10 frames (contact), else 6 A; the MHC residue sits across the periodic boundary
    box = np.array([50.0, 50.0, 50.0, 90, 90, 90])
    tg = [np.array([0, 1])] + [None] * 83
    mg = [np.array([2, 3, 4])] + [None] * 81
    f = PairMinDist(tg, mg)
    hits = 0
    for fr in range(10):
        X = np.array([[1.0, 0, 0], [1.0, 3, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
        gap = 4.0 if fr < 3 else 6.0
        X[2] = [1.0 - gap + 50.0, 0, 0]            # wrapped by +box: min image distance = gap
        X[3] = [1.0 - gap - 3 + 50.0, 0, 0]
        X[4] = [1.0, 0, 20.0]
        d = f(X, box)
        assert abs(d[0, 0] - gap) < 1e-5 and np.isnan(d[1, 1])   # MDAnalysis computes in float32
        hits += f.contacts(X, box)[0, 0]
    assert hits / 10 == 0.3
    # boundary: exactly 4.5 A is a contact (18b uses <=)
    X = np.array([[0.0, 0, 0], [0, 10, 0], [4.5, 0, 0], [30, 0, 0], [0, 30, 0]])
    assert f.contacts(X)[0, 0]


if __name__ == "__main__":
    tests = [(n, f) for n, f in sorted(globals().items()) if n.startswith("test_") and callable(f)]
    failed = 0
    for n, f in tests:
        try:
            f(); print(f"PASS  {n}")
        except AssertionError as e:
            failed += 1; print(f"FAIL  {n}  {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
