#!/usr/bin/env python3
"""Step 2 / gate G2: the contact function (contacts.PairMinDist) applied to each clean crystal PDB
must reproduce that complex's p82v1 label matrix Y exactly on resolved pairs (S = 1). Every mismatch
is reported; the cutoff and atom selection are NOT tuned.

Also checked here (supporting, not the gate):
  * brute force: an 18b-style per-pair loop gives the same minimum distances (to 1e-4 A);
  * PBC: wrapping the MHC chain by one box vector leaves minimum-image distances unchanged;
  * the same function on the MD starting coordinates (CHARMM-GUI step3_input, MD atom indices from
    the residue map) - informative only: CHARMM-GUI completes truncated side chains.
Reads: clean crystal PDBs, step3_input, residue maps, pair_maps_p82v1.npz (Y, S). No trajectory, no score.
usage: python md/analysis/validate_crystal_contacts.py
"""
import sys
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from mdmap import ROOT, SYSTEMS, crystal_model, chain_roles, crystal_slots, md_protein_universe  # noqa: E402
from contacts import PairMinDist, CUTOFF, groups_from_map  # noqa: E402

lab = np.load(ROOT / "data/processed/pair_maps_p82v1.npz", allow_pickle=True)
LPDB = list(lab["pdb"])
assert float(lab["cutoff"]) == CUTOFF


def crystal_groups(slots, coords):
    g = []
    for r in slots:
        if r is None:
            g.append(None); continue
        idx = []
        for a in r:
            if a.element != "H":
                idx.append(len(coords)); coords.append(a.coord.astype(np.float64))
        g.append(np.array(idx))
    return g


def brute(tcr, mhc):
    """18b's own loop, returning distances instead of a boolean."""
    out = np.full((84, 82), np.nan)
    for i, rt in enumerate(tcr):
        if rt is None: continue
        for j, rm in enumerate(mhc):
            if rm is None: continue
            out[i, j] = min(np.linalg.norm(a.coord - b.coord) for a in rt if a.element != "H"
                            for b in rm if b.element != "H")
    return out


all_ok = True
for pdb in SYSTEMS:
    i = LPDB.index(pdb)
    Y, S = lab["Y"][i].astype(bool), lab["S"][i].astype(bool)
    model = crystal_model(pdb)
    tcr, mhc = crystal_slots(model, chain_roles(pdb))
    coords = []
    tg, mg = crystal_groups(tcr, coords), crystal_groups(mhc, coords)
    X = np.array(coords)
    f = PairMinDist(tg, mg)
    d = f(X)
    Yh = np.nan_to_num(d, nan=np.inf) <= CUTOFF
    mism = np.argwhere((Yh != Y) & S)
    outside = int((Yh & ~S).sum())
    ok = len(mism) == 0 and outside == 0
    all_ok &= ok
    tl = lab["tcr_position_labels"]; mp = lab["mhc_positions"]

    # supporting checks
    ref = brute(tcr, mhc)
    brute_dev = float(np.nanmax(np.abs(ref - d)))
    mhc_atoms = np.concatenate([g for g in mg if g is not None])
    box = np.array([200.0, 200.0, 200.0, 90, 90, 90])
    Xw = X.copy(); Xw[mhc_atoms] += [200.0, 0.0, 0.0]       # MHC wrapped to the next image
    pbc_dev = float(np.nanmax(np.abs(f(Xw, box) - d)))
    near = np.sort(np.abs(d[S] - CUTOFF))[:3]

    print(f"\n=== {pdb}  G2 {'PASS' if ok else 'FAIL'} ===")
    print(f"  resolved pairs {int(S.sum())}; p82v1 contacts {int(Y[S].sum())}; recomputed {int(Yh[S].sum())}; "
          f"mismatches {len(mism)}; contacts outside S {outside}")
    for a, b in mism:
        print(f"    MISMATCH {tl[a]} x {mp[b]}: label {int(Y[a, b])}, min distance {d[a, b]:.3f} A")
    print(f"  brute-force (18b loop) max |diff| {brute_dev:.2e} A; PBC wrap max |diff| {pbc_dev:.2e} A; "
          f"closest distances to the cutoff: {', '.join(f'{x:.3f}' for x in near)} A from 4.5")

    # informative: MD starting coordinates with MD atom indices
    df = pd.read_csv(HERE / f"residue_map_{pdb}.tsv", sep="\t")
    u = md_protein_universe(pdb)
    fm = PairMinDist(*groups_from_map(df))
    dm = fm(u.atoms.positions.astype(np.float64))
    Ym = np.nan_to_num(dm, nan=np.inf) <= CUTOFF
    diff = np.argwhere((Ym != Y) & S)
    print(f"  [info] MD starting structure (step3_input, CHARMM-GUI-completed side chains): "
          f"{int(Ym[S].sum())} contacts, {len(diff)} differ from the label"
          + (": " + "; ".join(f"{tl[a]} x {mp[b]} label {int(Y[a, b])} -> {dm[a, b]:.2f} A" for a, b in diff) if len(diff) else ""))

print(f"\nG2 overall: {'PASS' if all_ok else 'FAIL'}")
sys.exit(0 if all_ok else 1)
