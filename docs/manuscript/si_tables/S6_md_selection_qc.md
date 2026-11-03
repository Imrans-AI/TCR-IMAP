# Table S6. Molecular dynamics: system selection and per-replicate quality control

**Table S6.** (A) The four simulated TCR–pMHC complexes: the predeclared selection basis, structure and their standing
under the corrected (p84v2) analysis. (B) Quality control of each of the twelve 200 ns production replicates.

## A. Systems

| PDB | Stratum | Selected as | Δ (p82v1) | TCR-IMAP (p82v1) | Prior (p82v1) | Resolution (Å) | clu70 | Allele | Peptide | p84v2 rank | Δ (p84v2) | TCR-IMAP (p84v2) | Prior (p84v2) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5yxn | high-Δ | pre-identified anchor | +0.2503 | 0.5360 | 0.2856 | 2.03 | 3 | HLA-A*02:01 | KLVALGINAV | 12 / 174 | +0.2446 | 0.5250 | 0.2804 |
| 3vxs | high-Δ | next-ranked eligible | +0.2167 | 0.6196 | 0.4029 | 1.80 | 28 | HLA-A*24:01 | RYPLTLGWCF | 4 / 174 | +0.2662 | 0.6660 | 0.3998 |
| 7n1f | low-Δ | pre-identified anchor | −0.2314 | 0.4490 | 0.6804 | 2.39 | 42 | not annotated | YLQPRTFLL | 167 / 174 | −0.1105 | 0.5694 | 0.6799 |
| 8eo8 | low-Δ | next-ranked eligible | −0.1946 | 0.1695 | 0.3641 | 2.30 | 38 | HLA-B*35:01 | LPFDKATIM | 154 / 174 | −0.0651 | 0.2853 | 0.3504 |

## B. Production replicates

| PDB | Rep | Steps | Frames (20 ps) | LINCS warnings | mdrun sessions | T (K), mean ± SD | P (bar), mean ± SD | Density (kg m⁻³), mean (drift 20–200 ns) | Backbone RMSD mean (Å) | Backbone RMSD max (Å) | TCR–pMHC contacts t = 0 → 200 ns (retained) | Docked % (count / identity) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5yxn | 1 | 50,000,000 | 10,001 | 0 | 3 | 300.07 ± 0.51 | 1.04 ± 41.8 | 1027.37 (+0.04) | 4.57 | 7.12 | 305 → 64 (21.0 %) | 5 / 0 |
| 5yxn | 2 | 50,000,000 | 10,001 | 0 | 3 | 300.06 ± 0.51 | 1.07 ± 42.0 | 1027.40 (+0.07) | 3.62 | 7.35 | 297 → 88 (29.6 %) | 52 / 0 |
| 5yxn | 3 | 50,000,000 | 10,001 | 0 | 2 | 300.06 ± 0.51 | 1.00 ± 41.8 | 1027.39 (−0.02) | 4.55 | 7.38 | 317 → 109 (34.4 %) | 3 / 0 |
| 7n1f | 1 | 50,000,000 | 10,001 | 0 | 2 | 300.06 ± 0.52 | 1.19 ± 42.2 | 1027.44 (−0.02) | 4.83 | 7.67 | 297 → 171 (57.6 %) | 78 / 0 |
| 7n1f | 2 | 50,000,000 | 10,001 | 0 | 3 | 300.06 ± 0.52 | 1.10 ± 42.2 | 1027.43 (+0.10) | 3.47 | 5.63 | 275 → 159 (57.8 %) | 93 / 11 |
| 7n1f | 3 | 50,000,000 | 10,001 | 0 | 1 | 300.06 ± 0.51 | 0.91 ± 42.3 | 1027.41 (−0.00) | 3.63 | 5.59 | 303 → 225 (74.3 %) | 100 / 22 |
| 3vxs | 1 | 50,000,000 | 10,001 | 0 | 1 | 300.06 ± 0.49 | 0.98 ± 39.8 | 1026.55 (−0.03) | 2.84 | 4.65 | 326 → 230 (70.6 %) | 100 / 12 |
| 3vxs | 2 | 50,000,000 | 10,001 | 0 | 2 | 300.06 ± 0.49 | 0.93 ± 39.9 | 1026.54 (+0.08) | 2.95 | 4.93 | 330 → 281 (85.2 %) | 100 / 19 |
| 3vxs | 3 | 50,000,000 | 10,001 | 0 | 1 | 300.06 ± 0.49 | 0.98 ± 39.9 | 1026.56 (−0.07) | 3.56 | 5.32 | 303 → 248 (81.8 %) | 100 / 29 |
| 8eo8 | 1 | 50,000,000 | 10,001 | 0 | 1 | 300.07 ± 0.53 | 1.05 ± 43.3 | 1028.01 (+0.02) | 3.60 | 5.58 | 285 → 209 (73.3 %) | 100 / 2 |
| 8eo8 | 2 | 50,000,000 | 10,001 | 0 | 1 | 300.07 ± 0.53 | 0.97 ± 43.4 | 1028.02 (+0.10) | 3.64 | 6.08 | 279 → 165 (59.1 %) | 98 / 0 |
| 8eo8 | 3 | 50,000,000 | 10,001 | 0 | 1 | 300.07 ± 0.53 | 0.94 ± 43.4 | 1028.00 (+0.03) | 4.15 | 6.48 | 271 → 203 (74.9 %) | 99 / 0 |

**Notes.**
- Δ = per-complex pairwise AUPRC of TCR-IMAP − pair-frequency prior under tcr70. **Selection used the p82v1 values**
  (the analysis version current when the criteria were committed), mechanically, before any system was built or
  simulated: the two highest- and two lowest-Δ complexes passing the predeclared filters (resolution ≤ 2.5 Å; distinct
  clu70 clusters; complete α1/α2 window; a single peptide of standard residues, no covalent TCR–peptide link, engineered
  interchain disulfide or single-chain construct; member of the 174-complex set). 5yxn and 7n1f were pre-identified
  anchors, included because they passed every filter. The selection was not revisited after the p84v2 correction; the
  p84v2 rank (1 = largest Δ) is given for transparency. The p84v2 extremes are 3qfj (+0.279) and 3qeq (−0.158).
  7n1f's p84v2 Δ is −0.110 (stored −0.110461). The criteria document's revision first gave −0.111, a double rounding
  through −0.1105, and was corrected to −0.110 on 28 Oct 2026; the rank, 167, is unchanged.
- Every replicate: 50,000,000 steps of 4 fs (200 ns), checkpoint complete, 10,001 frames at 20 ps with no irregular
  timesteps, 0 LINCS warnings. Interruptions were resumed from checkpoint, which the protocol allows. Temperature,
  pressure and density are the QC script's averages over the production run, with the density drift over 20–200 ns.
- Backbone RMSD: protein Cα RMSD from the energy-minimized structure over the analyzed window (20–200 ns, every frame,
  9,001 frames), on trajectories made whole and jump-free (the characterization recomputed after the periodic-boundary
  fix). The QC script's own coarser RMSD (every 10 ns, nm) is kept in the CSV.
- TCR–pMHC contacts: heavy-atom pairs < 4.5 Å between TCRα/TCRβ and the peptide plus MHC heavy chain (β2m excluded),
  summed over the four interfaces, at t = 0 and 200 ns; retained = 200 ns / t = 0, computed here from the QC log and
  given to one decimal. The criteria document's retention table gives the same values: its replicate-2 entries,
  recorded as integers in QC round 2 (17 Oct), were corrected to these log values on 28 Oct 2026.
- Docked %: share of the 9,001 window frames counted as docked, under the count reading (at least half as many TCR–pMHC
  heavy-atom contacts as at t = 0; the predeclared rule) and the identity reading (at least half of the specific t = 0
  contact pairs present).
- mdrun sessions are as logged by the QC script. They differ from the narrative in the criteria document ("All 12
  replicates against the predeclared protocol"), which names interruptions only for 5yxn rep1, 5yxn rep3, 7n1f rep2 and
  3vxs rep2; the log also records 3 sessions for 5yxn rep2 and 2 for 7n1f rep1. No other quantity is affected.

**Sources.** A: `md/MD_SYSTEM_SELECTION_CRITERIA.md` §8 (selection record: p82v1 Δ, TCR-IMAP, prior, resolution, clu70,
allele, peptide) and "Revision — the analysis moves to the p84v2 pair universe" (p84v2 ranks); p84v2 rank and values
recomputed from `results/tables/percomplex_tcr70_paired_p84v2.tsv` and asserted equal to that revision; clu70 and
resolution checked against `data/processed/complex_metadata.tsv`. B: `results/review/characterisation_fix/S2_qc_all12.log`
(`md/qc_replicate.sh` on all 12 replicates: steps, frames, LINCS, sessions, thermodynamics, interface contacts; retention
checked against the criteria document's retention table), `results/tables/md_characterisation_summary.csv` and
`md_characterisation_timeseries.csv` (backbone RMSD), `results/tables/md_occupancy_docked_frames.csv` (docked %).
