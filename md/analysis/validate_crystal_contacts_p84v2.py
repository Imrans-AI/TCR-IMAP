#!/usr/bin/env python3
"""Gate G2 on the p84v2 axes: the contact function (contacts.PairMinDist, unchanged) applied to each clean
crystal PDB must reproduce that complex's **p84v2** label matrix Y exactly on resolved pairs (S = 1). Every
mismatch is reported; the cutoff and atom selection are NOT tuned.

Also checked here (supporting, not the gate):
  * brute force: an 18b-style per-pair loop gives the same minimum distances;
  * PBC: wrapping the MHC chain by one box vector leaves minimum-image distances unchanged;
  * the same function on the MD starting coordinates (step3_input, MD atom indices from the p84v2 residue
    map) - informative only: CHARMM-GUI completes truncated side chains.

Reads: clean crystal PDBs, step3_input, residue_map_<pdb>_p84v2.tsv, pair_maps_p84v2.npz (Y, S).
No trajectory, no score.
usage: python md/analysis/validate_crystal_contacts_p84v2.py
"""
import sys
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from mdmap import ROOT, SYSTEMS, crystal_model, chain_roles, md_protein_universe  # noqa: E402
from mdmap_p84v2 import (crystal_slots_p84v2, groups_from_map_p84v2, LPDB, N_ROW, N_COL,  # noqa: E402
                         ROW_LABEL, COL_LABEL, _Z)
from contacts import PairMinDist, CUTOFF  # noqa: E402

assert float(_Z["cutoff"]) == CUTOFF


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
    out = np.full((N_ROW, N_COL), np.nan)
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
    Y, S = _Z["Y"][i].astype(bool), _Z["S"][i].astype(bool)
    model = crystal_model(pdb)
    tcr, mhc = crystal_slots_p84v2(pdb, model, chain_roles(pdb))
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

    ref = brute(tcr, mhc)
    brute_dev = float(np.nanmax(np.abs(ref - d)))
    mhc_atoms = np.concatenate([g for g in mg if g is not None])
    box = np.array([200.0, 200.0, 200.0, 90, 90, 90])
    Xw = X.copy(); Xw[mhc_atoms] += [200.0, 0.0, 0.0]
    pbc_dev = float(np.nanmax(np.abs(f(Xw, box) - d)))
    near = np.sort(np.abs(d[S] - CUTOFF))[:3]

    print(f"\n=== {pdb}  G2 (p84v2) {'PASS' if ok else 'FAIL'} ===")
    print(f"  resolved pairs {int(S.sum())}; p84v2 contacts {int(Y[S].sum())}; recomputed {int(Yh[S].sum())}; "
          f"mismatches {len(mism)}; contacts outside S {outside}")
    for a, b in mism:
        print(f"    MISMATCH {ROW_LABEL[a]} x {COL_LABEL[b]}: label {int(Y[a, b])}, min distance {d[a, b]:.3f} A")
    print(f"  brute-force (18b loop) max |diff| {brute_dev:.2e} A; PBC wrap max |diff| {pbc_dev:.2e} A; "
          f"closest distances to the cutoff: {', '.join(f'{x:.3f}' for x in near)} A from 4.5")
    ins = [j for j, lab in enumerate(COL_LABEL) if not lab.isdigit()]
    print(f"  insertion columns: " + "; ".join(
        f"{COL_LABEL[j]} {int(Y[:, j].sum())} contacts" for j in ins)
        + "; base columns " + "; ".join(f"{COL_LABEL[j]} {int(Y[:, j].sum())} contacts"
                                        for j in (COL_LABEL.index("1061"), COL_LABEL.index("1072"))))

    df = pd.read_csv(HERE / f"residue_map_{pdb}_p84v2.tsv", sep="\t")
    u = md_protein_universe(pdb)
    fm = PairMinDist(*groups_from_map_p84v2(df))
    dm = fm(u.atoms.positions.astype(np.float64))
    Ym = np.nan_to_num(dm, nan=np.inf) <= CUTOFF
    diff = np.argwhere((Ym != Y) & S)
    print(f"  [info] MD starting structure (step3_input, CHARMM-GUI-completed side chains): "
          f"{int(Ym[S].sum())} contacts, {len(diff)} differ from the label"
          + (": " + "; ".join(f"{ROW_LABEL[a]} x {COL_LABEL[b]} label {int(Y[a, b])} -> {dm[a, b]:.2f} A"
                              for a, b in diff) if len(diff) else ""))

print(f"\nG2 (p84v2) overall: {'PASS' if all_ok else 'FAIL'}")
sys.exit(0 if all_ok else 1)
