#!/usr/bin/env python3
"""Gate G1 on the p84v2 axes: write md/analysis/residue_map_<pdb>_p84v2.tsv and chain_segments_p84v2.tsv
for the four MD systems and check that the mapping is exact.

G1 passes for a system iff
  (a) every aligned crystal/MD residue pair has the same residue type,
  (b) every row/column residue present in the crystal is mapped to an MD residue (none in a gap),
  (c) each mapped MD residue sits on its crystal residue after one rigid superposition (CA within
      COORD_TOL), which resolves ambiguity between identical neighbours that sequence alone cannot -
      in particular Ala1061/Ala1061A, whose insertion codes CHARMM-GUI drops,
  (d) the coverage reproduces the p84v2 resolved-pair mask S of that complex exactly.

The p82v1 maps (residue_map_<pdb>.tsv) are NOT overwritten; these are separate files.
Reads: clean crystal PDBs, CHARMM-GUI step3_input.psf/.pdb, pair_maps_p84v2.npz (S and the stored
residue ids). No trajectory, no score.
usage: python md/analysis/build_residue_maps_p84v2.py
"""
import sys
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mdmap import ROOT, SYSTEMS, COORD_TOL  # noqa: E402
from mdmap_p84v2 import build_map_p84v2, LPDB, N_ROW, N_COL, _Z  # noqa: E402

OUT = ROOT / "md/analysis"
all_ok, chain_rows = True, []
for pdb in SYSTEMS:
    df, report, segids, roles, mdres, fit_rmsd = build_map_p84v2(pdb)
    df.to_csv(OUT / f"residue_map_{pdb}_p84v2.tsv", sep="\t", index=False)
    chain_rows += report
    S = _Z["S"][LPDB.index(pdb)]
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

    print(f"\n=== {pdb}  G1 (p84v2) {'PASS' if ok else 'FAIL'} ===")
    print("  chains: " + ", ".join(f"{k} {v}->{segids[v]}" for k, v in roles.items()))
    for r in report:
        print(f"  {r['role']:4s} {r['chain']}->{r['segid']}: crystal {r['n_crystal']}, MD {r['n_md']}, "
              f"aligned {r['n_aligned']}, type mismatches {r['mismatches']}, "
              f"crystal-only [{r['crystal_only'] or '-'}], MD-only (modelled) [{r['md_only'] or '-'}]")
    print(f"  superposition of all aligned CA: rmsd {fit_rmsd:.4f} A; mapped residues: CA deviation max "
          f"{mapped.ca_dev.max():.4f} A, heavy-atom deviation max {mapped.max_heavy_dev.max():.4f} A")
    print(f"  TCR rows with a crystal residue {int((t.status != 'absent_in_crystal').sum())}/{N_ROW}, "
          f"mapped {int((t.status == 'mapped').sum())}; MHC columns "
          f"{int((m.status != 'absent_in_crystal').sum())}/{N_COL}, mapped {int((m.status == 'mapped').sum())}")
    print(f"  resolved pairs: rebuilt S {int(S_rebuilt.sum())} == p84v2 S {int(S.sum())}: {s_equal}; "
          f"every resolved pair mapped: {s_mapped_equal}")
    ins = df[(df.status != "absent_in_crystal") & (df.crystal_icode.fillna("") != "")]
    print(f"  insertion-code residues on their OWN row/column: {len(ins)}")
    for _, r in ins.iterrows():
        print(f"    {r.kind} {r.slot_label} -> {r.crystal_chain} {r.crystal_resname}{int(r.crystal_resnum)}"
              f"{r.crystal_icode} -> MD {r.md_segid} {r.md_resname} resid {int(r.md_resid)}, residue "
              f"#{int(r.md_ordinal)} of the segment, CA dev {r.ca_dev:.4f} A")
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

pd.DataFrame(chain_rows).to_csv(OUT / "chain_segments_p84v2.tsv", sep="\t", index=False)
print(f"\nwrote md/analysis/residue_map_<pdb>_p84v2.tsv x {len(SYSTEMS)}, md/analysis/chain_segments_p84v2.tsv")
print(f"G1 (p84v2) overall: {'PASS' if all_ok else 'FAIL'}")
sys.exit(0 if all_ok else 1)
