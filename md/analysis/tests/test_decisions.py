#!/usr/bin/env python3
"""Synthetic tests with known answers for everything added by the 24 October 2026 decisions:
the identity-based docked reading, the distance-stratified rho, and the characterisation measures
(SASA, RMSD, RMSF). No MD trajectory and no model score is read.

usage: python md/analysis/tests/test_decisions.py
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import characterise as ch          # noqa: E402
import stats as st                 # noqa: E402
from occupancy import is_docked, is_docked_identity, DOCKED_FRAC  # noqa: E402

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' - ' + detail) if detail else ''}")


print("=== SASA (Shrake-Rupley, probe 1.4 A) ===")
from scipy.spatial import cKDTree  # noqa: E402
sph = ch.sphere_points()
X1, R1 = np.zeros((1, 3)), np.array([1.70])
s1 = ch._sasa_of(X1, R1, np.array([0]), cKDTree(X1), sph)[0]
check("isolated atom SASA = 4*pi*(r+probe)^2", abs(s1 - 4 * np.pi * (1.70 + 1.4) ** 2) < 1e-9,
      f"{s1:.4f} vs {4 * np.pi * (1.70 + 1.4) ** 2:.4f}")

# two identical atoms far apart: each keeps its full sphere
X2 = np.array([[0., 0, 0], [100., 0, 0]]); R2 = np.array([1.70, 1.70])
s2 = ch._sasa_of(X2, R2, np.arange(2), cKDTree(X2), sph)
check("two distant atoms are each fully exposed", np.allclose(s2, s1, atol=1e-9))

# two overlapping equal spheres, centres d apart: the buried part of each is a spherical cap of height
# h = R - d/2 (R = r + probe), area 2*pi*R*h, so the accessible area is 4*pi*R^2 - 2*pi*R*h = pi*R*(2R + d).
# (The fully coincident case d = 0 is degenerate - every point then lies exactly ON the neighbour's
# surface - and Shrake-Rupley's usual convention counts a point on the surface as accessible, so it is
# not a meaningful test of the implementation.)
r_, d_ = 1.70, 3.10
Rr = r_ + ch.SASA_PROBE
X3 = np.array([[0., 0, 0], [d_, 0, 0]]); R3 = np.array([r_, r_])
s3 = ch._sasa_of(X3, R3, np.arange(2), cKDTree(X3), sph)
want = np.pi * Rr * (2 * Rr + d_)
check("two overlapping spheres match the analytic cap formula (256 points, ~0.5 % quantisation)",
      np.allclose(s3, want, rtol=1e-2), f"{s3[0]:.3f} vs {want:.3f}")
s3k = ch._sasa_of(X3, R3, np.arange(2), cKDTree(X3), ch.sphere_points(4096))
# The error falls with refinement but not monotonically: the cap boundary cuts through the point
# lattice, so the quadrature error oscillates. 4096 points reach ~0.1 %.
check("and the error falls as the sphere is refined (4096 points, ~0.1 %)",
      np.allclose(s3k, want, rtol=2e-3) and abs(s3k[0] - want) < abs(s3[0] - want),
      f"{s3k[0]:.4f} vs {want:.4f}")

# buried SASA: groups far apart bury nothing
Xa = np.random.default_rng(0).normal(0, 3, (40, 3)); Xb = Xa + np.array([200.0, 0, 0])
X = np.vstack([Xa, Xb]); R = np.full(80, 1.70)
b_far = ch.buried_sasa(X, R, np.arange(40), np.arange(40, 80), sph)
check("buried SASA of two distant groups = 0", abs(b_far) < 1e-9, f"{b_far:.3e}")

# buried SASA: the interface-restricted computation equals the definition
rng = np.random.default_rng(1)
Xa = rng.normal(0, 4, (60, 3)); Xb = rng.normal(0, 4, (60, 3)) + np.array([6.0, 0, 0])
X = np.vstack([Xa, Xb]); R = np.full(120, 1.70)
fast = ch.buried_sasa(X, R, np.arange(60), np.arange(60, 120), sph)
naive = ch.buried_sasa_naive(X, R, np.arange(60), np.arange(60, 120), sph)
check("buried SASA (interface-restricted) = buried SASA (definition)", abs(fast - naive) < 1e-6,
      f"{fast:.4f} vs {naive:.4f}")
check("buried SASA of touching groups is positive", fast > 1.0, f"{fast:.1f} A^2")

print("\n=== RMSD / RMSF ===")
A = rng.normal(0, 5, (50, 3))
check("RMSD of a structure with itself = 0", ch.rmsd(A, A) == 0.0)
th = 0.7
Rz = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
B = A @ Rz.T + np.array([10.0, -3, 2])
check("fitted RMSD removes a rigid-body move", ch.fitted_rmsd(B, A, np.arange(50), np.arange(50)) < 1e-9)
# RMSF of a known oscillation: x alternates +-d about the mean -> RMSF = d
d = 0.4
acc = ch.RMSFAccumulator(3)
for k in range(200):
    P = np.zeros((3, 3)); P[:, 0] = d * (1 if k % 2 == 0 else -1)
    acc.add(P)
check("RMSF of a +-d square wave = d", np.allclose(acc.rmsf(), d, atol=1e-12), f"{acc.rmsf()}")

print("\n=== docked-frame readings (decision 3) ===")
check("count-based: exactly 50 % of the t = 0 count is docked", is_docked(50, 100))
check("count-based: 49 % is not", not is_docked(49, 100))
check("identity-based: exactly 50 % of the t = 0 pairs is docked", is_docked_identity(50, 100))
check("identity-based: 49 % is not", not is_docked_identity(49, 100))
# the two readings can disagree: same count, completely different pairs
ref = np.arange(100)
now = np.arange(1000, 1100)          # 100 contacts again, none of them the original ones
kept = int(np.isin(ref, now).sum())
check("a frame can satisfy the count rule and fail the identity rule",
      is_docked(len(now), len(ref)) and not is_docked_identity(kept, len(ref)),
      f"count {len(now)}/{len(ref)}, persisting {kept}")
check("DOCKED_FRAC is the predeclared 0.50", DOCKED_FRAC == 0.50)

print("\n=== distance-stratified rho (decision 4) ===")
n = 60
Y = np.zeros((n, 1)); S = np.ones((n, 1))
dist = np.zeros((n, 1)); score = np.zeros((n, 1)); occ = np.zeros((n, 1))
# stratum 1 (4.5-6 A): perfectly increasing; stratum 2 (6-8): perfectly decreasing; stratum 3 (>8): flat-ish
for k in range(20):
    dist[k, 0] = 5.0; score[k, 0] = k; occ[k, 0] = k
for k in range(20, 40):
    dist[k, 0] = 7.0; score[k, 0] = k; occ[k, 0] = -k
for k in range(40, 60):
    dist[k, 0] = 12.0; score[k, 0] = k; occ[k, 0] = (k % 3)
rows = {r["stratum"]: r for r in st.rho_by_stratum(score, occ, Y, S, dist, n_perm=200)}
check("stratum 4.5-6 A: rho = +1", abs(rows["4.5-6 A"]["rho"] - 1.0) < 1e-12, f"{rows['4.5-6 A']['rho']}")
check("stratum 6-8 A: rho = -1", abs(rows["6-8 A"]["rho"] + 1.0) < 1e-12, f"{rows['6-8 A']['rho']}")
check("each stratum has its own n = 20", all(rows[k]["n"] == 20 for k in rows))
# boundaries: 6.0 belongs to 6-8, not to 4.5-6
db = np.array([[4.5], [5.9999], [6.0], [7.9999], [8.0], [100.0]])
Yb = np.zeros((6, 1)); Sb = np.ones((6, 1))
masks = st.stratum_masks(Yb, Sb, db)
check("4.5 A is in the near stratum", masks["4.5-6 A"][0, 0])
check("6.0 A falls in 6-8, not 4.5-6", masks["6-8 A"][2, 0] and not masks["4.5-6 A"][2, 0])
check("8.0 A falls in > 8", masks["> 8 A"][4, 0] and not masks["6-8 A"][4, 0])
check("every crystal-negative pair lands in exactly one stratum",
      np.array_equal(sum(masks.values()).astype(int), np.ones((6, 1), int)))
# a small stratum is reported with its n and no test
small = st.rho_by_stratum(score[:10], occ[:10], Y[:10], S[:10], dist[:10], n_perm=50)
check("a stratum with fewer than 20 pairs is reported untested",
      np.isnan(small[0]["rho"]) and small[0]["n"] == 10 and "not tested" in small[0]["note"])

print(f"\n{len(PASS)}/{len(PASS) + len(FAIL)} synthetic checks passed")
if FAIL:
    print("FAILED: " + ", ".join(FAIL))
sys.exit(1 if FAIL else 0)
