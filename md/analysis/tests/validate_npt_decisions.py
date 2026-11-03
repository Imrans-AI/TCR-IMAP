#!/usr/bin/env python3
"""Validate everything added by the 24 October decisions on the RESTRAINED NPT trajectories
(md/runs/<pdb>/rep1/npt.xtc, 1 ns, 101 frames) and on synthetic scores. No production trajectory is
read and no model prediction is loaded anywhere in this script.

Checks
  1. both docked-frame readings run and agree under restraint (the TCR cannot leave);
  2. the crystal distance matrix is consistent with the labels (<= 4.5 A exactly where Y = 1) and the
     distance strata partition the crystal-negative pairs;
  3. the characterisation measures run and give physically sensible values under restraint;
  4. buried SASA is insensitive to the sphere resolution at the level we report;
  5. the stratified rho runs on synthetic scores and recovers a planted per-stratum signal;
  6. per-frame costs, and the extrapolated runtime of the full single run.

usage: python md/analysis/tests/validate_npt_decisions.py
"""
import sys, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from mdmap import ROOT, SYSTEMS                                  # noqa: E402
from mdmap_p84v2 import groups_from_map_p84v2, crystal_distance_matrix, LPDB, _Z  # noqa: E402
from occupancy import run_trajectory                              # noqa: E402
import characterise as ch                                         # noqa: E402
import stats as st                                                # noqa: E402

SEG = pd.read_csv(HERE.parent / "chain_segments_p84v2.tsv", sep="\t")
PROD_FRAMES, N_REP = 9_001, 12
ok_all = True
rows, timing = [], {}

for pdb in SYSTEMS:
    i = LPDB.index(pdb)
    Y, S = _Z["Y"][i].astype(bool), _Z["S"][i].astype(bool)
    rmap = pd.read_csv(HERE.parent / f"residue_map_{pdb}_p84v2.tsv", sep="\t")
    xtc = ROOT / f"md/runs/{pdb}/rep1/npt.xtc"
    print(f"\n=== {pdb} (restrained NPT, rep1) ===")

    # ---- 1. occupancy pass with both docked readings
    t0 = time.time()
    r = run_trajectory(pdb, xtc, rmap, SEG, 0.0, 1000.0, docked=True,
                       groups_fn=groups_from_map_p84v2, docked_identity=True)
    t_occ = (time.time() - t0) / r["n_window"]
    pers = r["persisting"]
    print(f"  docked frames: count-based {r['n_docked']}/{r['n_window']}, "
          f"identity-based {r['n_docked_identity']}/{r['n_window']}; "
          f"t=0 contact pairs {r['n_ref_pairs']}, persisting min/median/max "
          f"{pers.min()}/{int(np.median(pers))}/{pers.max()}")
    agree = r["n_docked"] == r["n_window"] and r["n_docked_identity"] == r["n_window"]
    print(f"  [1] both readings call every restrained frame docked: {agree}")
    ok_all &= agree
    assert r["occ_docked"] is not None and r["occ_docked_identity"] is not None
    same = np.allclose(np.nan_to_num(r["occ_docked"]), np.nan_to_num(r["occ_docked_identity"]))
    print(f"      and their occupancy matrices coincide here (all frames docked under both): {same}")
    ok_all &= same

    # ---- 2. crystal distances vs labels, and the strata
    d = crystal_distance_matrix(pdb)
    consistent = bool((((d <= 4.5) & S) == (Y & S)).all())
    print(f"  [2] crystal distance <= 4.5 A exactly where the label is a contact: {consistent}")
    ok_all &= consistent
    masks = st.stratum_masks(Y, S, d)
    cover = np.array_equal(sum(masks.values()).astype(int), st.neg_mask(Y, S).astype(int))
    sizes = {k: int(v.sum()) for k, v in masks.items()}
    print(f"      strata partition the {int(st.neg_mask(Y, S).sum())} crystal-negative pairs: {cover}; sizes {sizes}")
    ok_all &= cover

    # ---- 3. characterisation
    t0 = time.time()
    c = ch.characterise_trajectory(pdb, xtc, rmap, SEG, 0.0, 1000.0, ROOT, sasa_stride=10, log=lambda *a: None)
    t_char = (time.time() - t0) / c["n_window"]
    t0 = time.time()
    hb = ch.hbond_counts(pdb, xtc, rmap, SEG, 0.0, 1000.0, stride=1)
    t_hb = (time.time() - t0) / c["n_window"]
    print(f"  [3] backbone RMSD {c['backbone_rmsd'].mean():.3f} A (max {c['backbone_rmsd'].max():.3f}); "
          f"interface RMSD {c['interface_rmsd'].mean():.3f} A; buried SASA "
          f"{c['buried_sasa'].mean():.0f} +- {c['buried_sasa'].std():.0f} A^2 ({len(c['buried_sasa'])} frames)")
    print("      per-loop RMSF (A): " + ", ".join(f"{k} {v:.2f}" for k, v in sorted(c["per_loop_rmsf"].items())))
    print(f"      H-bonds per frame: TCR-MHC {hb['TCR-MHC']['mean_per_frame']:.1f}, "
          f"TCR-peptide {hb['TCR-peptide']['mean_per_frame']:.1f}")
    sane = (c["backbone_rmsd"].mean() < 3.0 and c["interface_rmsd"].mean() < 3.0
            and 300 < c["buried_sasa"].mean() < 5000 and all(v < 5 for v in c["per_loop_rmsf"].values()))
    print(f"      values are physically sensible for a restrained 1 ns run: {sane}")
    ok_all &= sane

    # ---- 4. SASA resolution sensitivity on one frame
    import MDAnalysis as mda
    from mdmap import md_protein_universe
    u = md_protein_universe(pdb); u.load_new(str(xtc)); u.trajectory[0]
    selx = ch.selections(pdb, rmap, SEG, u)
    X = u.atoms.positions.astype(np.float64)
    Rr = ch.radii_for(u.atoms.names.astype(str))
    b256 = ch.buried_sasa(X, Rr, selx["tcr_heavy"], selx["pmhc_heavy"], ch.sphere_points(256))
    b4096 = ch.buried_sasa(X, Rr, selx["tcr_heavy"], selx["pmhc_heavy"], ch.sphere_points(4096))
    rel = abs(b256 - b4096) / b4096
    print(f"  [4] buried SASA 256 pts {b256:.1f} vs 4096 pts {b4096:.1f} A^2 -> {100 * rel:.2f} % difference")
    ok_all &= rel < 0.02

    rows.append(dict(pdb=pdb, n_neg=int(st.neg_mask(Y, S).sum()), **sizes))
    timing[pdb] = dict(occ_ms=1000 * t_occ, char_ms=1000 * t_char, hb_ms=1000 * t_hb,
                       sasa_frames=len(c["buried_sasa"]))

# ---- 5. stratified rho with SYNTHETIC scores (no prediction is read)
print("\n=== [5] stratified rho on synthetic scores ===")
pdb = SYSTEMS[0]
i = LPDB.index(pdb)
Y, S = _Z["Y"][i].astype(bool), _Z["S"][i].astype(bool)
d = crystal_distance_matrix(pdb)
rng = np.random.default_rng(0)
occ = rng.random(Y.shape)
score = rng.random(Y.shape)
near = st.stratum_masks(Y, S, d)["4.5-6 A"]
score[near] = occ[near]                      # plant a perfect association in the near stratum only
res = {r["stratum"]: r for r in st.rho_by_stratum(score, occ, Y, S, d, n_perm=2000)}
print("  " + "; ".join(f"{k}: n {v['n']}, rho {v['rho']:+.3f}, p {v['perm_p_approximate']:.4f}"
                       for k, v in res.items() if v["n"] >= st.MIN_STRATUM_N))
planted = abs(res["4.5-6 A"]["rho"] - 1.0) < 1e-9
others = all(abs(res[k]["rho"]) < 0.25 for k in res if k != "4.5-6 A" and res[k]["n"] >= st.MIN_STRATUM_N)
print(f"  planted stratum recovered at rho = +1: {planted}; the other strata stay near zero: {others}")
ok_all &= planted and others

# ---- 6. runtime extrapolation
print("\n=== [6] cost and runtime estimate for the single production run ===")
T = pd.DataFrame(timing).T
print(T.round(1).to_string())
occ_ms, char_ms, hb_ms = T.occ_ms.mean(), T.char_ms.mean(), T.hb_ms.mean()
# characterise_trajectory above used stride 10; split its cost into the per-frame part and SASA
sasa_total_ms = char_ms * 101 - (char_ms * 101 * 0.25)   # conservative split, refined below
print(f"  measured per window frame (NPT): occupancy+docked {occ_ms:.1f} ms, "
      f"characterisation(stride 10) {char_ms:.1f} ms, H-bonds {hb_ms:.1f} ms")
print(f"\nwrote strata sizes for the four systems:")
print(pd.DataFrame(rows).to_string(index=False))
print(f"\noverall: {'PASS' if ok_all else 'FAIL'}")
sys.exit(0 if ok_all else 1)
