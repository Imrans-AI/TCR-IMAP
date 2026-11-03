#!/usr/bin/env python3
"""Gate G3 on the p84v2 axes: run the occupancy code on md/runs/<pdb>/rep1/npt.xtc (restrained NPT, 1 ns,
10 ps frames) for all four systems. Equilibration data only; no score is read.

G3 expectation (task): crystal-positive pairs show high occupancy under the 400/40 kJ/mol/nm2
restraints. Reported: the distribution of crystal-positive occupancy over all 101 frames, the
crystal-negative pairs that reach any contact (mechanics check), the TCR-pMHC contact count
relative to t = 0 (exercises the docked-frame code), frame timing, and per-frame cost.
G3 is judged PASS if every system's crystal-positive median occupancy is >= 0.9; the full
distribution is printed whatever the outcome.
usage: python md/analysis/validate_npt.py
"""
import sys
import time
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from mdmap import ROOT, SYSTEMS  # noqa: E402
from occupancy import run_trajectory
from mdmap_p84v2 import groups_from_map_p84v2  # noqa: E402

lab = np.load(ROOT / "data/processed/pair_maps_p84v2.npz", allow_pickle=True)
LPDB = list(lab["pdb"]); TL = lab["row_label"]; MP = lab["col_label"]
seg = pd.read_csv(HERE / "chain_segments_p84v2.tsv", sep="\t")

all_ok, rows = True, []
for pdb in SYSTEMS:
    i = LPDB.index(pdb)
    Y, S = lab["Y"][i].astype(bool), lab["S"][i].astype(bool)
    rmap = pd.read_csv(HERE / f"residue_map_{pdb}_p84v2.tsv", sep="\t")
    xtc = ROOT / f"md/runs/{pdb}/rep1/npt.xtc"
    t0 = time.time()
    r = run_trajectory(pdb, xtc, rmap, seg, t_start_ps=0.0, t_end_ps=1000.0, docked=True, groups_fn=groups_from_map_p84v2)
    dt = time.time() - t0
    occ = r["occ"]
    pos = occ[Y & S]; neg = occ[~Y & S]
    assert not np.isnan(occ[S]).any(), "NaN occupancy on a resolved pair"
    step = np.diff(r["times"])
    ok = np.median(pos) >= 0.9
    all_ok &= ok
    print(f"\n=== {pdb} rep1 NPT  G3 (p84v2) {'PASS' if ok else 'FAIL'} ===")
    print(f"  frames {r['n_window']} (t {r['times'][0]:.0f}-{r['times'][-1]:.0f} ps, spacing "
          f"{step.min():.0f}-{step.max():.0f} ps); {dt:.1f} s total, {1000 * dt / r['n_window']:.1f} ms/frame")
    q = np.quantile(pos, [0, 0.1, 0.25, 0.5, 0.75, 1])
    print(f"  crystal-positive pairs n = {len(pos)}: occupancy min {q[0]:.2f}, p10 {q[1]:.2f}, q1 {q[2]:.2f}, "
          f"median {q[3]:.2f}, q3 {q[4]:.2f}, max {q[5]:.2f}; mean {pos.mean():.3f}; "
          f">= 0.25: {int((pos >= 0.25).sum())}, >= 0.5: {int((pos >= 0.5).sum())}, >= 0.9: {int((pos >= 0.9).sum())}")
    lo = np.argwhere((occ < 0.5) & Y & S)
    for a, b in lo:
        print(f"    low crystal-positive: {TL[a]} x {MP[b]} occupancy {occ[a, b]:.2f}")
    print(f"  crystal-negative pairs n = {len(neg)}: any contact {int((neg > 0).sum())}, >= 0.25 {int((neg >= 0.25).sum())}, "
          f">= 0.5 {int((neg >= 0.5).sum())}  (mechanics check under restraint; no scores involved)")
    c = r["counts"]
    print(f"  TCR-pMHC heavy-atom contacts: t = 0 {r['ref_count']}, window min {c.min()} / median {int(np.median(c))} / "
          f"max {c.max()}; docked frames (>= 50 % of t = 0) {r['n_docked']}/{r['n_window']}")
    rows.append(dict(pdb=pdb, n_frames=r["n_window"], n_pos=len(pos), pos_median=np.median(pos),
                     pos_mean=pos.mean(), pos_min=pos.min(), n_pos_ge_0_9=int((pos >= 0.9).sum()),
                     n_neg=len(neg), neg_any=int((neg > 0).sum()), neg_ge_0_25=int((neg >= 0.25).sum()),
                     ref_count=r["ref_count"], docked=r["n_docked"], ms_per_frame=1000 * dt / r["n_window"]))

print(f"\nG3 (p84v2) overall: {'PASS' if all_ok else 'FAIL'}")
print(pd.DataFrame(rows).round(3).to_string(index=False))
sys.exit(0 if all_ok else 1)
