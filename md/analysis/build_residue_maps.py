#!/usr/bin/env python3
"""Step 1 / gate G1: write md/analysis/residue_map_<pdb>.tsv and chain_segments.tsv for the four MD
systems and check that the mapping is exact.

G1 passes for a system iff
  (a) every aligned crystal/MD residue pair has the same residue type,
  (b) every slot residue present in the crystal is mapped to an MD residue (none in a gap),
  (c) each mapped MD residue sits on its crystal residue after one rigid superposition (CA within
      COORD_TOL; step3_input carries crystal coordinates) - this catches alignment ambiguity between
      identical neighbours, which sequence alone cannot,
  (d) the slot coverage reproduces the p82v1 resolved-pair mask S of that complex exactly.
Gaps (crystal-only or MD-only residues) are reported in full whether or not they touch a slot.

Reads: clean crystal PDBs, CHARMM-GUI step3_input.psf/.pdb, data/processed/pair_maps_p82v1.npz
(S only). No trajectory, no score.
usage: python md/analysis/build_residue_maps.py
"""
import sys
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mdmap import ROOT, SYSTEMS, COORD_TOL, build_map  # noqa: E402

OUT = ROOT / "md/analysis"
lab = np.load(ROOT / "data/processed/pair_maps_p82v1.npz", allow_pickle=True)
LPDB = list(lab["pdb"])

all_ok, chain_rows = True, []
for pdb in SYSTEMS:
    df, report, segids, roles, mdres, fit_rmsd = build_map(pdb)
    df.to_csv(OUT / f"residue_map_{pdb}.tsv", sep="\t", index=False)
    for r in report:
        chain_rows.append(r)
    S = lab["S"][LPDB.index(pdb)]
    t = df[df.kind == "tcr"].set_index("slot")
    m = df[df.kind == "mhc"].set_index("slot")
    S_rebuilt = np.outer(t.status != "absent_in_crystal", m.status != "absent_in_crystal").astype(int)
    S_mapped = np.outer(t.status == "mapped", m.status == "mapped").astype(int)

    mism = sum(r["mismatches"] for r in report)
    bad = df[df.status.isin(["NOT_IN_MD", "RESNAME_MISMATCH"])]
    mapped = df[df.status == "mapped"]
    coord_ok = bool((mapped.ca_dev <= COORD_TOL).all())
    s_equal = bool((S_rebuilt == S).all())
    s_mapped_equal = bool((S_mapped == S).all())
    ok = mism == 0 and len(bad) == 0 and coord_ok and s_equal and s_mapped_equal
    all_ok &= ok

    print(f"\n=== {pdb}  G1 {'PASS' if ok else 'FAIL'} ===")
    print("  chains: " + ", ".join(f"{k} {v}->{segids[v]}" for k, v in roles.items()))
    for r in report:
        print(f"  {r['role']:4s} {r['chain']}->{r['segid']}: crystal {r['n_crystal']}, MD {r['n_md']}, "
              f"aligned {r['n_aligned']}, type mismatches {r['mismatches']}, "
              f"crystal-only [{r['crystal_only'] or '-'}], MD-only (modelled) [{r['md_only'] or '-'}]")
    print(f"  superposition of all aligned CA: rmsd {fit_rmsd:.4f} A; mapped slot residues: CA deviation max "
          f"{mapped.ca_dev.max():.4f} A, heavy-atom deviation max {mapped.max_heavy_dev.max():.4f} A")
    print(f"  TCR slots with a crystal residue {int((t.status != 'absent_in_crystal').sum())}/84, "
          f"mapped {int((t.status == 'mapped').sum())}; MHC positions "
          f"{int((m.status != 'absent_in_crystal').sum())}/82, mapped {int((m.status == 'mapped').sum())}")
    print(f"  resolved pairs: rebuilt S {int(S_rebuilt.sum())} == p82v1 S {int(S.sum())}: {s_equal}; "
          f"every resolved pair mapped: {s_mapped_equal}")
    for _, r in df[(df.status != "absent_in_crystal") & (df.crystal_icode.fillna("") != "")].iterrows():
        print(f"  insertion-code residue on a slot: {r.kind} {r.slot_label} -> {r.crystal_chain} "
              f"{r.crystal_resname}{int(r.crystal_resnum)}{r.crystal_icode} (18b last-wins) -> MD "
              f"{r.md_segid} {r.md_resname} resid {int(r.md_resid)}, residue #{int(r.md_ordinal)} of the segment, "
              f"CA dev {r.ca_dev:.4f} A")
    for _, r in df[df.note.fillna("") != ""].iterrows():
        print(f"  alignment corrected by coordinates: {r.kind} {r.slot_label}: {r.note}")
    for _, r in bad.iterrows():
        print(f"  PROBLEM {r.status}: {r.kind} {r.slot_label} {r.crystal_chain} {r.crystal_resname}"
              f"{int(r.crystal_resnum)}{r.crystal_icode} {r.get('note', '')}")
    neq = mapped[mapped.n_heavy_crystal != mapped.n_heavy_md]
    print(f"  heavy-atom counts crystal = MD on {len(mapped) - len(neq)}/{len(mapped)} mapped residues"
          + ("; differ: " + ", ".join(f"{r.slot_label} {r.crystal_resname} {int(r.n_heavy_crystal)}->{int(r.n_heavy_md)}"
                                      for _, r in neq.iterrows()) if len(neq) else ""))
    shared = [x for x in mdres[segids[roles["MHC"]]] if sum(y.resid == x.resid for y in mdres[x.segid]) > 1]
    if shared:
        print("  MD residues sharing a resid (insertion codes dropped by CHARMM-GUI; split by atom-name repeat): "
              + ", ".join(f"{x.resname}{x.resid}" for x in shared))

seg = pd.DataFrame(chain_rows)
seg.to_csv(OUT / "chain_segments.tsv", sep="\t", index=False)
print(f"\nwrote md/analysis/residue_map_<pdb>.tsv x {len(SYSTEMS)}, md/analysis/chain_segments.tsv")
print(f"G1 overall: {'PASS' if all_ok else 'FAIL'}")
sys.exit(0 if all_ok else 1)
