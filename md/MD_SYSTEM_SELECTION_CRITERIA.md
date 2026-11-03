# MD system selection — predeclared criteria

**Written and committed BEFORE `percomplex_tcr70_paired_p82v1.tsv` was opened for selection.**
Git history is the evidence: this file's commit precedes any commit containing a selected PDB list.

Author: Muhammad Imran. Date: 2 October 2026.

---

## 1. Why this document exists

The primary analysis of this paper was rebuilt once already because a 28-position MHC candidate
set had been chosen using contact frequency pooled over all 174 complexes *before* splitting —
outcome-informed selection that invalidated every headline number (see `CLAUDE.md`, "The
correction"). The same failure mode is available here: choosing MD systems by looking at which
complexes give a convenient answer. These criteria are fixed in advance so that the selection
cannot be, and cannot be alleged to be, informed by the outcome.

## 2. Question the simulations address

Among residue pairs that TCR-IMAP scores highly but that are **negative in the crystal
structure**, what fraction show meaningful contact occupancy in solution?

If a real fraction of the model's apparent false positives are transient solution-phase contacts,
then the crystal-derived label set systematically penalises the model and measured performance is
a floor rather than an estimate. This addresses the paper's weakest point — the modest advantage
over the within-fold frequency prior — from an independent direction.

**Key statistic:** Spearman rho between TCR-IMAP score and MD contact occupancy, computed among
crystal-negative pairs only.

## 3. Predeclared outcome handling

**If high-scoring crystal-negative pairs show LOW MD occupancy, that result is reported.** It
would mean the model's false positives are genuinely wrong, and the paper says so. This is agreed
before any system is built so the analysis cannot be judged on its result. No re-running with
different systems, seeds, force fields or cutoffs in search of a different answer.

## 4. Selection criteria — fixed

Source table: `results/tables/percomplex_tcr70_paired_p82v1.tsv` (tcr70 protocol, per-complex
paired differences, p82v1).

Selection quantity: **Delta = per-complex pairwise AUPRC(TCR-IMAP) - AUPRC(pair-frequency prior)**
under tcr70. This is the comparison the MD analysis is designed to inform.

**Four systems, two strata:**
- **2 high-Delta** — complexes where the model most exceeds the prior
- **2 low-Delta** — complexes where it does not

The mid-Delta stratum in the earlier 6-system sketch (`CLAUDE.md`, MD design table) is dropped.
Four systems preserve two complexes per stratum, so no stratum rests on a single structure.

**Eligibility filters, applied before ranking:**
1. Crystallographic resolution <= 2.5 A
2. Distinct clu70 cluster — no two selected complexes from the same cluster
3. Complete alpha1/alpha2 platform: no unresolved residues inside the 82-position window
4. Single peptide, standard residues only; no covalent TCR-peptide linkage, no engineered
   interchain disulfide, no single-chain construct
5. Member of the 174-complex primary set (`6uz1` is out of scope by construction)

**Ranking rule:** among eligible complexes, take the 2 highest Delta and the 2 lowest Delta.

**Pre-identified extremes.** `5yxn` (Delta = +0.250) and `7n1f` (Delta = -0.231) are already
recorded in `CLAUDE.md` and `HANDOFF_1.md` as the Figure 6 success/failure extremes, so the ends
of this table are already known to the author. They are therefore **auto-included as the two
anchors if and only if they pass the eligibility filters above**; the second complex in each
stratum is the next-ranked eligible complex. If either fails a filter, it is replaced by the
next-ranked eligible complex and the failure is recorded here. This rule is stated now so that
prior knowledge of the extremes cannot be exercised as a choice later.

**Ties** are broken by higher resolution (lower numeric value), then by lexicographic PDB ID.

## 5. Simulation protocol — fixed

| Element | Setting |
|---|---|
| Builder | CHARMM-GUI **Solution Builder** (soluble complex; not Membrane Builder) |
| Force field | CHARMM36m, TIP3P water |
| Box | cubic, 10 A buffer, 150 mM KCl, neutralised |
| Engine | GROMACS 2023.2, CUDA (RTX 5070 Ti, driver 570.211.01) |
| Electrostatics | PME |
| Constraints | LINCS on H-bonds, **4 fs timestep with hydrogen mass repartitioning** |
| Equilibration | staged minimisation -> NVT 300 K -> NPT 300 K / 1 bar |
| **Replicates** | **3 per system, different initial velocity seeds** |
| **Production** | **200 ns per replicate, FIXED IN ADVANCE** (revised 2 Oct 2026, see below) |
| Output | protein-only XTC every 20 ps (~0.7 GB per replicate) |
| PBC | existing three-stage pipeline: whole -> nojump -> center |

Total: 4 systems x 3 replicates x 200 ns = **2,400 ns**.

### Protocol revision, 2 October 2026 — before any production run, before any analysis

Two changes, both decided from measured throughput and applied identically to all twelve
replicates. Nothing had been simulated for analysis and no occupancy number had been computed
when these were made.

**1. Hydrogen mass repartitioning, 4 fs timestep.** Hydrogen masses scaled x3 with the surplus
subtracted from the bonded heavy atom; water and ions untouched (SETTLE handles rigid TIP3P).
Hydrogens 1.008 -> 3.024 amu; mass conserved per molecule to < 1e-9 amu; lightest resulting
heavy atom 5.963 amu (a methyl carbon, 12.011 - 3 x 2.016), identical in all four systems.
*(Corrected 6 Oct 2026: an earlier version of this line gave 3.024 amu as the lightest heavy
atom. That was a repartitioned hydrogen, miscounted because the check treated any atom above
1.5 amu as heavy. No topology was affected; only the reported number was wrong.)*
The committed script reproduces the benchmarked `3vxs` topology with zero per-atom mass
difference, and was applied unchanged to all four systems. Build provenance (CHARMM-GUI job
IDs, bundle checksums, atom counts) is in `md/systems/BUILD_MANIFEST.tsv`. Measured on
`3vxs` (415,551 atoms): **77.3 ns/day at 2 fs -> 145.7 ns/day at 4 fs**, a 1.89x speedup.
Stability over a 2 ns test: mean temperature 300.07 K, RMSD 0.63 K, drift -0.23 K.
HMR preserves equilibrium structural and thermodynamic properties and alters only
hydrogen-dependent dynamics, which this analysis does not measure (Hopkins et al., *JCTC* 2015,
11, 1864). Script: `scripts/md_apply_hmr.py`; original topologies retained unmodified.

**2. Production length 500 ns -> 200 ns per replicate.** At the measured 145.7 ns/day, the
original 4 x 3 x 500 ns design required ~41 days of continuous GPU time. 200 ns per replicate
gives 2,400 ns total, ~17 days, and retains the full four-system two-stratum design rather than
buying time by cutting systems. Contact occupancy is a local property that converges well within
this window, and 200 ns exceeds the 3 x 100 ns replicate standard used for comparable analyses in
this journal.

**The design was NOT reduced in system count.** Dropping to two systems would have collapsed each
stratum to n = 1, making any stratum-level difference indistinguishable from between-structure
variation. Length was cut instead because it costs the least evidentially.

**Any future extension** (for example in response to peer review) applies to all twelve replicates
and is recorded here with its date and external trigger. Extension prompted by inspection of the
occupancy result is excluded by the predeclaration below.

**No trajectory is extended after its result is seen.** The 200 ns length is fixed here. If the
schedule forces a change, it applies to all twelve replicates and is recorded in this file with
its date and reason, before the occupancy analysis is run.

**All three replicates are pooled** for the occupancy statistic; per-replicate values are reported
in SI so the spread is visible. No "representative trajectory" is selected post hoc. (This departs
deliberately from the lab's TvTPPIMI protocol, which selected one stable replicate for extended
analysis — defensible for a ligand-stability narrative, but it would bias an occupancy estimate
upward by construction.) Replicate design follows Knapp, Ospina & Deane, *JCTC* 2018, 14,
6127-6138.

### Run protocol — fixed 6 October 2026, before launch

Files: `md/protocol/{nvt,npt,prod}.mdp`, `md/run_replicate.sh`, `md/run_queue.sh`. All four
CHARMM-GUI builds generated byte-identical mdp files, so one shared protocol applies to all.
The HMR topology (`topol_hmr.top`) is used in every stage.

| Stage | Length | dt | Ensemble | Restraints (BB/SC, kJ/mol/nm2) |
|---|---|---|---|---|
| Minimisation | CHARMM-GUI step4.0 | - | - | as generated; run once per system, shared by its replicates |
| NVT | 125 ps | 1 fs | v-rescale 300 K | 400 / 40; **velocities generated here, seeded** |
| NPT | 1 ns | 2 fs | + C-rescale 1 bar | 400 / 40; `refcoord-scaling = com` |
| Production | 200 ns | 4 fs | v-rescale 300 K, C-rescale 1 bar | none |

The restrained NPT stage is not part of CHARMM-GUI's pipeline; it is added so the run matches the
committed sequence (minimisation -> NVT -> NPT) and box density relaxes under restraint rather
than during production.

**Velocity seeds, fixed:** replicate 1 = 271828, replicate 2 = 314159, replicate 3 = 161803,
applied at the NVT stage of every system. These seeds are what make the three replicates
independent; they are recorded here so any replicate can be regenerated exactly.

**Output:** protein-only (`SOLU`) compressed trajectory every 20 ps -> 10,001 frames per replicate.
The first 20 ns of production is discarded before analysis (section 6).

**Execution order:** one replicate at a time, alternating strata -
`5yxn` (high), `7n1f` (low), `3vxs` (high), `8eo8` (low) - completing replicate 1 for all four
before starting replicate 2. An interruption therefore never leaves one stratum covered and the
other not. Every stage is resume-safe; production restarts from its checkpoint.

**Validation:** `gmx grompp -maxwarn 0` passed for every stage of every system (16 of 16) with
zero warnings. The same `-maxwarn 0` is used in the real runs, so no warning can be suppressed.

**Execution (revised 6 Oct 2026, before any replicate ran):** all twelve equilibrations run first;
production then runs in checkpointed sessions (overnight and weekends) and is resumed from
checkpoint, which continues the same trajectory without altering its ensemble. A stage counts as
complete only when its checkpointed step reaches `nsteps`; an interrupted `mdrun` still writes a
final structure, so output files are not treated as evidence of completion.

**Known integration note, accepted:** at 4 fs, grompp notes that side-chain carbonyl bonds
(Asp/Glu/Asn/Gln CG-OD1, CD-OE1; period ~22 fs) are sampled at fewer than 10 steps per period.
Neither atom carries hydrogen, so HMR does not affect them. This is inherent to 4 fs protocols
with hydrogen-only constraints (as validated by Hopkins et al. 2015) and is below the
velocity-Verlet stability limit (period / pi ~ 7 fs); the 2 ns benchmark held 300.07 K with
-0.23 K drift. Constraining all bonds instead would depart from how CHARMM36m was parameterised.

### Production QC — round 1 (13 October 2026), recorded before any outcome analysis

`md/qc_replicate.sh` checks completion, trajectory continuity, thermodynamics, backbone RMSD and
gross interface integrity. It deliberately computes no per-pair occupancy and nothing involving
TCR-IMAP scores (section 6 is run only once all 12 replicates exist).

All four replicate-1 runs are mechanically clean: 50,000,000 steps, 10,001 frames at 20 ps with no
irregular timesteps (5yxn rep1 crashed once at 177 ns and resumed seamlessly from its checkpoint),
0 LINCS warnings, T 300.06-300.07 K, P ~1 bar, density drift <= 0.04 kg/m3 over 20-200 ns. Peptide-MHC,
MHC-beta2m and TCRalpha-TCRbeta interfaces persist in every system.

**Observation:** share of starting TCR-pMHC heavy-atom contacts (< 4.5 A) still present at 200 ns -
8eo8 73 %, 3vxs 71 %, 7n1f 58 %, **5yxn 21 %**. In 5yxn rep1 the TCR substantially disengaged
(beta-chain-peptide 17 -> 0, beta-chain-MHC 94 -> 14) without fully dissociating.

**No action taken in response.** No replicate re-run, replaced or excluded, and no protocol change.
Whether the docked state is maintained is exactly what section 6 check 1 (occupancy of
crystal-positive pairs) measures; its result will be reported per replicate, and interpretation of
crystal-negative occupancy will account for it. The remaining 5yxn replicates show whether this is
characteristic of the system or a single-replicate excursion.

### Secondary analysis — docked frames (added 13 October 2026, before any occupancy was computed)

Added after the round-1 integrity QC showed the TCR disengaging in 5yxn rep1. It is therefore
informed by that observation and is reported as a **secondary sensitivity analysis only**. The
primary analysis (section 6: all frames after the 20 ns discard, all replicates pooled) is unchanged
and remains the headline result.

- **Docked frame:** a frame retaining >= 50 % of the TCR-pMHC heavy-atom contacts (< 4.5 A) present in
  the first production frame (t = 0) of the same replicate. TCR-pMHC = TCRalpha and TCRbeta against
  the peptide and the MHC heavy chain (beta2m excluded), summed - the definition used by
  `md/qc_replicate.sh`. The t = 0 reference is the end of restrained NPT, i.e. the crystal docking mode.
- **Statistic:** identical to the primary - Spearman rho between TCR-IMAP score and MD occupancy among
  crystal-negative pairs - computed on docked frames only, per system and pooled.
- **Minimum data:** a replicate contributes only if it has >= 20 ns (1,000 frames) of docked frames.
- **One threshold, fixed at 50 %.** No other threshold will be examined or reported.
- **Also reported:** the fraction of docked frames for every replicate, whatever the outcome.

At the time this was written, no per-pair occupancy and no quantity involving TCR-IMAP scores had
been computed from any trajectory.

### Production QC — round 2 (17 October 2026), recorded before any outcome analysis

All four replicate-2 runs are mechanically clean: 50,000,000 steps, 10,001 frames at 20 ps with no
irregular timesteps, 0 LINCS warnings, T 300.06-300.07 K, P ~1 bar, density drift <= 0.10 kg/m3,
backbone RMSD plateaued at 0.29-0.44 nm. 7n1f rep2 (one segfault plus a machine reboot, 3 mdrun
sessions) and 3vxs rep2 (one shutdown, 2 sessions) show no trace of the interruptions in the data.

**TCR-pMHC contact retention at 200 ns (share of the t = 0 count):**

| System | rep1 | rep2 |
|---|---|---|
| 3vxs | 70.6 % | 85.2 % |
| 7n1f | 57.6 % | 57.8 % |
| 8eo8 | 73.3 % | 59.1 % |
| **5yxn** | **21.0 %** | **29.6 %** |

*(Corrected 28 Oct 2026: three replicate-2 values were recorded here as 57 % (7n1f), 58 % (8eo8) and 29 % (5yxn).
The QC log (`results/review/characterisation_fix/S2_qc_all12.log`, the four TCR<->peptide/MHC heavy-atom contact
counts summed, 200 ns / t = 0) gives 57.8 %, 59.1 % and 29.6 %. All values are now given to one decimal from that
log; the replicate-1 values agree with the integers recorded on 13 October. No conclusion changes.)*

5yxn is low in both replicates while the other three stay in the 57.6-85.2 % band, with the same
detail in each: the beta chain loses nearly all contacts (rep2 peptide 18 -> 2, MHC 86 -> 17) while
the alpha chain partially recovers by 200 ns. The round-1 observation is therefore not a
single-replicate excursion; rep3 will confirm whether it holds across all three.

**No action taken.** No replicate re-run, replaced or excluded; no protocol change. Section 6
check 1 reports per-replicate crystal-positive occupancy, and the docked-frame sensitivity analysis
predeclared on 13 October applies its fixed 50 % threshold and >= 20 ns minimum as written - under
which 5yxn may contribute few or no docked frames. That is the rule operating as intended.

## 5b. System construction — fixed before any system was built

**Chains retained:** MHC-I heavy chain (full ectodomain, alpha1/alpha2/alpha3), **beta-2
microglobulin**, peptide, TCR alpha and TCR beta variable+constant domains as resolved.

beta2m and the alpha3 domain are **retained**, not trimmed. They are structural components of the
native MHC-I heterotrimer; deleting them would expose large buried surfaces, destabilise alpha3,
and risk propagating artefacts into the alpha1/alpha2 platform that the analysis measures. The
cost is roughly 100 residues on a ~400-residue complex.

**Removed before building:** crystallographic waters, cryoprotectants, buffer ions and any
crystallisation additives; for `8eo8`, the additional copies in the asymmetric unit (20 chains
present) - the single biological assembly A/C/D/E is used.

**Missing residues** are modelled by CHARMM-GUI only where they fall inside a chain; disordered
terminal tails are not extended. Any complex requiring modelling of residues inside the
82-position window would be recorded here - none do, since all four have 82/82 resolved.

Protonation states at pH 7.4 as assigned by CHARMM-GUI; histidine tautomers left at the default
assignment and not hand-tuned per system.

## 6. Analysis — fixed

Contact occupancy uses the **same minimum heavy-atom 4.5 A definition as the p82v1 labels**, so
MD occupancy and crystal labels are directly comparable.

1. Occupancy of each **crystal-positive** pair (expected high; sanity check)
2. Occupancy of each **crystal-negative** pair, split by TCR-IMAP score
3. **Spearman rho between TCR-IMAP score and occupancy among crystal-negative pairs** — the key
   statistic, per system and pooled
4. Per-loop RMSF, germline (CDR1/2) vs CDR3
5. Interface RMSD; buried SASA

Occupancy is computed over the post-equilibration production window only; the discarded initial
window is fixed at 20 ns per replicate and is not tuned.

## 7. Stated limits

- Four systems is a **case-controlled analysis, not a benchmark**, and is described as such.
- MD occupancy is force-field dependent and is **not ground truth either**. It is a second
  independent view of the same interface, not an arbiter.
- Agreement between TCR-IMAP scores and MD occupancy would not establish that the model learns
  physical residue-residue interaction; that claim stays withdrawn.

## 8. Selection record

Table opened 2 October 2026, after criteria commit `268ceed`.

**Filter 3 (complete 82-position window) was applied and excluded nothing** - all 174
complexes have all 82 window positions resolved. Recorded for transparency; it was vacuous.

Filter 1 (resolution <= 2.5 A) left 75 of 174 complexes across 31 clu70 clusters. Filter 4 was
checked against the structure files: peptide chain length, residue standardness, chain
distinctness and absence of engineered interchain disulfides. Both pre-identified anchors
passed every filter and are auto-included per section 4. No replacements were required.

| Stratum | PDB | Delta | TCR-IMAP | Prior | Res | clu70 | Allele | Peptide | Selected as |
|---|---|---|---|---|---|---|---|---|---|
| high-Delta | **5yxn** | +0.2503 | 0.5360 | 0.2856 | 2.03 A | 3 | HLA-A*02:01 | KLVALGINAV | pre-identified anchor |
| high-Delta | **3vxs** | +0.2167 | 0.6196 | 0.4029 | 1.80 A | 28 | HLA-A*24:01 | RYPLTLGWCF | next-ranked eligible |
| low-Delta | **7n1f** | -0.2314 | 0.4490 | 0.6804 | 2.39 A | 42 | (not annotated) | YLQPRTFLL | pre-identified anchor |
| low-Delta | **8eo8** | -0.1946 | 0.1695 | 0.3641 | 2.30 A | 38 | HLA-B*35:01 | LPFDKATIM | next-ranked eligible |

All four occupy distinct clu70 clusters (3, 28, 42, 38), distinct alleles and distinct source
pathogens (HCV NS3, HIV Nef, SARS-CoV-2 spike, influenza NP). Peptides are two 10-mers and two
9-mers, all standard residues, no non-standard or engineered components.

`5jzi` (Delta +0.2101) ranked third at the high end but shares clu70 cluster 3 with `5yxn` and
is excluded by the distinct-cluster rule; `3vxs` outranks it regardless.

Note on `7n1f`: the pair-frequency prior attains 0.6804 here, its highest value anywhere in the
set. The model is not failing so much as being outperformed on an unusually predictable
interface - an informative low-Delta case rather than a degenerate one.

Selection was performed by applying the committed criteria mechanically. No complex was chosen,
substituted or excluded outside the stated rules, and the two-strata design of section 4 is
unchanged.

**Status:** systems selected, none built, nothing simulated.

## Analysis definitions - fixed 20 October 2026, before any production occupancy was computed

- **Frames:** production 20-200 ns, all 9,001 frames per replicate (20 ps spacing).
- **Occupancy of a pair:** fraction of those frames with minimum heavy-atom distance <= 4.5 A (the
  p82v1 label definition). System-level occupancy = mean of its three replicate occupancies.
- **Pair universe:** the 84 x 82 candidate pairs resolved in the crystal labels (S = 1).
  Crystal-positive: Y = 1. Crystal-negative: Y = 0.
- **Scores:** TCR-IMAP tcr70 out-of-fold pairwise scores, as used in Figures 4 and 6.
- **Primary statistic:** Spearman rho between score and occupancy among crystal-negative pairs, per
  system and pooled over all four systems' crystal-negative pairs; per-replicate rho reported
  alongside.
- **Uncertainty:** within-system permutation p (10,000 permutations of scores among crystal-negative
  pairs), explicitly labelled approximate because pairs are not independent.
- **Meaningful occupancy:** >= 0.25, fixed. **High-scoring crystal-negative pairs:** the top k by
  score, k = that complex's number of crystal-positive pairs (the Figure 6 rule). Reported: the
  fraction of these with occupancy >= 0.25, against the base rate among all crystal-negative pairs.
- **Control, added here before any occupancy exists:** the same statistics computed with the
  within-fold pair-frequency prior scores instead of TCR-IMAP. Without it, a positive rho for
  TCR-IMAP cannot be distinguished from general positional contact propensity.
- **Check 1 (section 6):** crystal-positive occupancy per replicate. Docked-frame secondary analysis
  exactly as already predeclared (13 October).

At the time this section was committed, no per-pair occupancy, contact frequency, distance statistic
or correlation had been computed from any production trajectory, and no analysis code existed.

## Revision — the analysis moves to the p84v2 pair universe and predictions, 24 October 2026

Recorded **before any production occupancy was computed**. At the time of writing, no per-pair occupancy,
contact frequency, distance statistic or correlation had been computed from any production trajectory,
`md/analysis/RUN_AUTHORISED` did not exist, and no `results/tables/md_occupancy_*.csv` had been written.

### What changes

The analysis defined on 20 October named the **p82v1** pair universe (84 TCR slots x 82 MHC positions,
residues keyed by IMGT number alone) and the tcr70 out-of-fold scores trained on it. Both are replaced by
their **p84v2** equivalents:

- **Pair universe:** the p84v2 label array, **174 x 75 x 84** — 75 TCR rows (68 occupied IMGT positions plus
  7 CDR3 insertion rows) and 84 MHC columns (82 base positions plus **1061A** and **1072A**), every residue
  keyed by (IMGT number, insertion code). Crystal-positive is Y = 1 and crystal-negative Y = 0 in that array,
  over its resolved pairs S = 1, exactly as before.
- **Scores:** the tcr70 out-of-fold TCR-IMAP predictions trained on the p84v2 inputs
  (`results/tables/pairmaps_tcr70_p84v2.npz`), with the within-fold pair-frequency prior from the same file
  as the control.
- **Reason:** the p82v1 representation mis-assigned residues. Label column "1061" held residue 1061A in all
  174 complexes while residue 1061 itself was never labelled, and packed TCR tokens were not aligned with
  their label rows. The correction and its verification gates are recorded in
  `results/review/CORRECTION_PLAN_P84V2.md` and `REBUILD_P84V2_REPORT.md`.

**Everything else in this document is unchanged**: the question, the predeclared outcome handling, the
occupancy definition (fraction of retained frames with minimum heavy-atom distance <= 4.5 A), the 20-200 ns
window, the primary statistic (Spearman rho among crystal-negative pairs, per system and pooled), the
permutation procedure, the >= 0.25 meaningful-occupancy threshold, the prior control, check 1, and the
docked-frame secondary analysis with its fixed 50 % threshold and >= 20 ns minimum.

### System selection is NOT revisited

The four systems were selected mechanically on the **p82v1** per-complex deltas, after the criteria were
committed and before any system was built. **That selection stands exactly as predeclared and is not
re-run.** Re-selecting on the corrected deltas after seeing them would be the outcome-informed selection
this document exists to prevent.

For transparency: under p84v2 the four systems rank 12th (`5yxn`, delta +0.245), 4th (`3vxs`, +0.266),
167th (`7n1f`, -0.110) and 154th (`8eo8`, -0.065) of 174 by the same quantity. They remain two clearly
high-delta and two clearly low-delta complexes, so the two-strata design is intact, but they are no longer
the mechanical extremes. The extremes under p84v2 are `3qfj` (+0.279) and `3qeq` (-0.158).
*(Corrected 28 Oct 2026: 7n1f was given as -0.111, a double rounding of the stored -0.110461
(`percomplex_tcr70_paired_p84v2.tsv`) through -0.1105. The rank is unaffected.)*

### Gates re-passed on the p84v2 axes (24 October 2026)

The residue map, contact function and occupancy mechanics were rebuilt for the new axes and re-validated.
No gate was adjusted. Full numbers: `results/review/FINISH_PENDING_REPORT.md`.

| Gate | Result |
|---|---|
| **G1** residue mapping | **PASS** in all four systems. 84/84 MHC columns mapped (was 82/82): **1061A and 1072A now map to their own columns**, and base 1061/1072 are mapped separately. TCR rows with a crystal residue 59/55/60/65 of 75, all mapped. Superposition RMSD 0.0005 A, maximum CA deviation 0.0008 A. Rebuilt S equals the p84v2 S exactly in all four |
| **G2** contact function vs the p84v2 crystal labels | **PASS**, 0 mismatches in all four (5yxn 36 contacts, 7n1f 22, 3vxs 30, 8eo8 31). The counts rise against p82v1 for 5yxn (35) and 8eo8 (29) because base 1061/1072 now have their own columns |
| **G3** occupancy on the restrained NPT trajectories | **PASS**. Crystal-positive median occupancy 1.00, 1.00, 1.00, 0.99; means 0.87-0.93; every NPT frame docked |

### Still open — decisions for the author, not taken here

These four were raised on 20 October and remain open. **No code change has decided any of them**; the
implemented behaviour is stated so the choice is explicit.

1. **Permutation sidedness.** Implemented: two-sided. The question in section 2 is directional, so a
   one-sided test is defensible. One argument in `pipeline.statistics_phase`.
2. **Which "top k" pairs.** Implemented literally: the top k *crystal-negative* pairs by score, k = that
   complex's number of crystal-positive pairs. The other reading is the crystal-negative members of the
   model's top-k set over all resolved pairs, which is fewer than k pairs.
3. **"Docked frame".** Implemented: the TCR-pMHC heavy-atom contact *count* is >= 50 % of the t = 0 count.
   The alternative reading is that >= 50 % of the *specific* t = 0 atom pairs persist, which is stricter.
4. **Per-loop RMSF, interface RMSD and buried SASA** (section 6, items 4-5) still have no fixed definitions
   and are not built. They need their own definitions before any production data is touched, or they are
   dropped from the paper.

### Status

All 12 replicates are complete (every `prod.cpt` at step 50,000,000). The analysis is built, re-validated on
the p84v2 axes and **has not been run**. The run script refuses unless `md/analysis/RUN_AUTHORISED` exists;
only the author creates that file.

## Analysis decisions settled, and the characterisation measures added — 24 October 2026

Recorded **before any production occupancy was computed** and before any code implementing them exists. At the time
of writing, no per-pair occupancy, contact frequency, distance statistic, correlation, RMSD, RMSF, SASA or hydrogen
bond had been computed from any production trajectory; `md/analysis/RUN_AUTHORISED` did not exist and no
`results/tables/md_occupancy_*.csv` had been written.

**Everything in this section is reported in the paper or its SI regardless of outcome.** Nothing here is a screen
to be passed before a result is shown.

### Primary statistic — unchanged

Spearman rho between the TCR-IMAP score and MD occupancy among **crystal-negative resolved pairs**, per system and
pooled over the four systems, with per-replicate values reported alongside. The pair universe, the occupancy
definition, the 20-200 ns window, the >= 0.25 meaningful-occupancy threshold and the pair-frequency-prior control
are as already fixed.

### The four open decisions, now settled

1. **The permutation test is TWO-SIDED.** 10,000 permutations of the scores among each system's crystal-negative
   pairs, within system for the pooled statistic, seed 42. It is labelled **approximate** wherever it is reported,
   because candidate pairs are not independent: they share residues, and neighbouring pairs share atoms. The
   question in section 2 is directional, but a two-sided test is the conservative choice and is not conditioned on
   the direction the data turn out to take.

2. **"Top k" keeps its original meaning:** the **top k crystal-negative pairs by score**, with k set to that
   complex's number of crystal-positive pairs - the Figure 6 rule, stable descending sort in row-major pair order,
   ties at the boundary flagged. Unchanged from the 20 October definition.

3. **"Docked frame" is reported under BOTH readings**, and both are reported whether or not they agree:
   - **(a) count-based, primary** - the frame retains at least 50 % of the *number* of TCR-pMHC heavy-atom contacts
     present at t = 0 of that replicate. This is the rule predeclared on 13 October and the quantity
     `md/qc_replicate.sh` computes.
   - **(b) identity-based** - at least 50 % of the *specific* heavy-atom contact pairs present at t = 0 are still
     present in the frame. This is the stricter reading of the same English sentence.
     *(Corrected the same day, before any code was written: an earlier version of this line said "residue-residue
     contact pairs". The count-based rule (a) counts heavy-atom pairs, as `md/qc_replicate.sh` does, and the open
     question this decision answers was phrased in terms of atom pairs, so the identity-based reading must track
     the same objects for the two to be comparable. No code or number is affected; nothing had been run.)*
   Both docked-frame fractions, and the docked-frame secondary rho under each, are reported. Where they disagree
   the disagreement is itself reported; neither is chosen after the fact.

4. **Distance-stratified rho is added (NEW).** *Why it is needed:* a crystal-negative pair whose residues are
   already close in the crystal structure can show high MD occupancy through proximity alone. A positive pooled rho
   could therefore reflect interface geometry - the model scoring near-interface pairs highly because they are
   near the interface - rather than genuine transient contact. To separate these, rho is reported **within strata of
   the crystal-frame minimum heavy-atom distance** of each crystal-negative pair:

   | Stratum | Crystal-frame minimum heavy-atom distance |
   |---|---|
   | near | 4.5 - 6 A (just outside the contact cutoff) |
   | middle | 6 - 8 A |
   | far | > 8 A |

   The stratum boundaries are fixed here and will not be adjusted. Strata with fewer than 20 crystal-negative pairs
   in a system are reported with their n and no p value. The pooled value is reported alongside the strata, and the
   **same stratification is applied to the pair-frequency-prior control**, so that any stratum-level difference
   between model and prior can be read the same way as the pooled one.

### Simulation quality and interface characterisation

All of the following are computed in the single authorised run and **reported in full**. They are **descriptive**.
With four systems they characterise what happened in these simulations; they are **not** used to claim any
relationship between loop flexibility, interface burial or hydrogen bonding and prediction quality, and no such
claim will be made from them.

- **Backbone RMSD** - protein Ca RMSD against the energy-minimised starting structure, per replicate, over the
  20-200 ns window.
- **Per-loop RMSF** - per-residue Ca RMSF over 20-200 ns after a least-squares fit on the MHC alpha1/alpha2 Ca
  atoms, averaged within each of the six CDR loops, reported per replicate and as a per-system mean.
- **Interface RMSD** - Ca RMSD of the TCR CDR-loop residues after fitting on the MHC alpha1/alpha2 Ca atoms,
  relative to the crystal frame, per replicate.
- **Buried SASA at the TCR-pMHC interface** - SASA(TCR) + SASA(pMHC) - SASA(complex) per frame; mean and standard
  deviation over the window, and the time course.
- **Interface hydrogen bonds** - donor-acceptor heavy-atom distance <= 3.5 A with D-H...A angle >= 150 degrees,
  counted between TCR and MHC alpha1/alpha2 and between TCR and peptide **separately**; mean count per frame and
  per-residue-pair hydrogen-bond occupancy.

**5yxn's TCR disengagement is a known feature of these simulations and is to be characterised by these measures,
not explained away.** Production QC recorded the share of t = 0 TCR-pMHC contacts still present at 200 ns as 21.0 %
for replicate 1 (round 1, 13 October) and 29.6 % for replicate 2 (round 2, 17 October; recorded then as 29 %,
corrected 28 Oct 2026 from the QC log), against 57.6-85.2 % for the other
three systems, with the beta chain losing nearly all contacts while the alpha chain partially recovers. *A round-3
QC record was never written into this document* - section "Production QC - round 2" says only that replicate 3
would show whether the behaviour holds across all three - so whether it does is not established in writing here.
Check 1 (crystal-positive occupancy per replicate), the docked-frame fractions under both readings of decision 3,
and the characterisation measures above will settle it in the single run, and will be reported for 5yxn exactly as
for the other three systems. No replicate is excluded, re-run or down-weighted on account of it.

### Production QC — round 3 (24 October 2026), recorded before any outcome analysis

Run with `md/qc_replicate.sh <pdb> 3` for all four systems; full output in
`results/review/pre_md_run/A_qc_round3.log`. As in rounds 1 and 2 this computes no per-pair occupancy and
nothing involving TCR-IMAP scores. Recorded late: round 3 finished earlier but was never written into this
document, so statements about replicate 3 previously rested on nothing durable.

| | 5yxn | 7n1f | 3vxs | 8eo8 |
|---|---|---|---|---|
| checkpointed step | 50,000,000 | 50,000,000 | 50,000,000 | 50,000,000 |
| frames (20 ps), irregular steps | 10,001, 0 | 10,001, 0 | 10,001, 0 | 10,001, 0 |
| LINCS warnings | 0 | 0 | 0 | 0 |
| mdrun sessions (queue events) | 2 (1 CRASH, resumed) | 1 | 1 | 1 |
| temperature (K) | 300.06 +- 0.51 | 300.06 +- 0.51 | 300.06 +- 0.49 | 300.07 +- 0.53 |
| pressure (bar) | 1.00 +- 41.8 | 0.91 +- 42.3 | 0.98 +- 39.9 | 0.94 +- 43.4 |
| density (kg/m3), drift 20-200 ns | 1027.4, -0.02 | 1027.4, -0.00 | 1026.6, -0.07 | 1028.0, +0.03 |
| backbone RMSD (nm): mean 100-200 ns / max / final | 0.44 / 0.74 / 0.41 | 0.40 / 0.45 / 0.40 | 0.34 / 0.46 / 0.31 | 0.46 / 0.54 / 0.40 |

**Interface heavy-atom contacts (< 4.5 A), start -> min -> end:**

| Interface | 5yxn | 7n1f | 3vxs | 8eo8 |
|---|---|---|---|---|
| pep-MHC | 309 -> 204 -> 220 | 319 -> 270 -> 273 | 411 -> 294 -> 357 | 297 -> 234 -> 254 |
| pep-TCRa | 79 -> 5 -> 14 | 96 -> 75 -> 99 | 71 -> 42 -> 63 | 6 -> 4 -> 7 |
| pep-TCRb | 24 -> 5 -> 7 | 79 -> 30 -> 34 | 53 -> 43 -> 72 | 77 -> 50 -> 75 |
| TCRa-MHC | 120 -> 0 -> 20 | 105 -> 55 -> 85 | 95 -> 30 -> 41 | 118 -> 44 -> 58 |
| TCRb-MHC | 94 -> 31 -> 68 | 23 -> 7 -> 7 | 84 -> 55 -> 72 | 70 -> 13 -> 63 |
| MHC-b2m | 422 -> 276 -> 344 | 412 -> 221 -> 288 | 435 -> 302 -> 368 | 443 -> 285 -> 306 |
| TCRa-TCRb | 683 -> 459 -> 565 | 623 -> 415 -> 451 | 641 -> 429 -> 436 | 692 -> 489 -> 631 |

### All 12 replicates against the predeclared protocol

**All 12 meet it**: every replicate reached the full 50,000,000 steps, wrote 10,001 frames at 20 ps with no
irregular timesteps, and logged 0 LINCS warnings. The only deviations from an uninterrupted run are session
counts, and the protocol explicitly allows resumption from checkpoint: 5yxn rep1 (crash at 177 ns), 5yxn rep3
(1 crash, 2 sessions), 7n1f rep2 (segfault plus a machine reboot, 3 sessions) and 3vxs rep2 (one shutdown,
2 sessions). None shows any trace of the interruption in its trajectory, thermodynamics or frame spacing.

### 5yxn across all three replicates — what the numbers show

Share of the t = 0 TCR-pMHC heavy-atom contacts (TCRalpha and TCRbeta against peptide and MHC heavy chain,
summed; b2m excluded) still present at 200 ns:

| System | rep1 | rep2 | rep3 |
|---|---|---|---|
| 3vxs | 70.6 % | 85.2 % | **81.8 %** |
| 7n1f | 57.6 % | 57.8 % | **74.3 %** |
| 8eo8 | 73.3 % | 59.1 % | **74.9 %** |
| **5yxn** | **21.0 %** | **29.6 %** | **34.4 %** |

*(Corrected 28 Oct 2026: the replicate-2 column carried the round-2 integers 57 % (7n1f), 58 % (8eo8) and 29 % (5yxn);
the QC log `results/review/characterisation_fix/S2_qc_all12.log` gives 57.8 %, 59.1 % and 29.6 %. All twelve values are
now given to one decimal from that log; the others agree with the integers previously recorded.)*

**5yxn loses most of its TCR-pMHC interface in all three replicates.** Its retention (21.0 %, 29.6 %, 34.4 %) is below
every replicate of every other system (57.6-85.2 %), and the three values are consistent with each other, so this is a
property of the system and not a single-replicate excursion. That answers the question left open in the round-2
record.

**But which chain disengages is NOT consistent across replicates**, and the round-2 wording should not be
generalised. In replicates 1 and 2 the beta chain lost nearly all of its contacts while the alpha chain partially
recovered. **In replicate 3 the opposite happens:** TCRalpha-MHC falls 120 -> 0 -> 20 and peptide-TCRalpha
79 -> 5 -> 14, while TCRbeta-MHC holds at 94 -> 31 -> 68. So the interface is substantially disrupted every time,
but the surviving chain differs between replicates. 5yxn rep3 also has the largest backbone RMSD excursion of the
twelve (max 0.74 nm, settling to 0.41 nm).

No replicate is excluded, re-run, replaced or down-weighted on account of any of this. Check 1 (crystal-positive
occupancy per replicate), the docked-frame fractions under both readings and the characterisation measures
predeclared on 24 October will quantify it in the single authorised run, and will be reported for 5yxn exactly as
for the other three systems.

## Undefined correlations are reported as undefined — 24 October 2026, before any production occupancy was seen

Recorded **before any implementing code was written**, and before any production occupancy has been written to disk,
displayed or summarised.

**The defect being corrected.** As implemented, when the occupancy values in a group are all identical, Spearman rho
is undefined (0/0). The permutation p was nevertheless returned as 1/(1 + 10,000) = 0.0001: every comparison against
an undefined observed value is false, so no permutation counts as "at least as extreme". A correlation that does not
exist was therefore reported with the smallest p the test can produce. The same happens for a group of fewer than two
pairs. The end-to-end test on the restrained NPT trajectories (`results/review/MD_CRASH_FIX_REPORT.md`) found it in
all ten "> 8 A" stratum rows, where every pair has occupancy exactly 0 under restraint. In production it is plausible
in any system where every far pair stays at zero occupancy.

**The rule.** A rank correlation is **undefined** for a group when:
- the group has fewer than two pairs, or
- the occupancy values in the group are all identical (including all zero), or
- the scores in the group are all identical.

Where it is undefined, **rho and p are both reported as undefined (NaN)**. No permutation test is run for that group.
**An undefined correlation is never reported as significant**, is never imputed as zero, and is not counted as
evidence in either direction.

So the reader can see why a row is undefined, every affected table reports the group's n and the number of distinct
occupancy values in the group (`n_distinct_occupancy`):

| Table | Change |
|---|---|
| `md_occupancy_rho.csv` | already carries n (`n_crystal_negative`); gains `n_distinct_occupancy` |
| `md_occupancy_strata.csv` | already carries `n`; gains `n_distinct_occupancy`, and its `note` states the reason for an undefined row |
| `md_occupancy_docked_rho.csv` | gains `n_crystal_negative` and `n_distinct_occupancy` |

**What does not change.** This changes only how an undefined correlation is reported. The following are all as
already fixed:
- the statistic (Spearman rho among crystal-negative resolved pairs);
- the pair universe, the occupancy definition and the 20–200 ns window;
- the permutation scheme: within system, 10,000 permutations, seed 42, two-sided, approximate;
- the distance strata and their boundaries, and the 20-pair minimum below which a stratum is reported untested;
- the >= 0.25 threshold, the top-k rule, the docked-frame readings and their 1,000-frame contribution rule;
- the prior control, and which rows each table contains.

A group of exactly two pairs with distinct values has a defined rho of +1 or −1; its two-sided permutation p is 1,
and it is reported as computed.

**Timing, stated exactly.** The first authorised attempt at the single run (authorised 24 October 2026, 09:12 UTC)
completed the 12 Phase A trajectory passes and then crashed in `occupancy_tables` on a container-type defect
(`results/review/MD_CRASH_FIX_REPORT.md`):
- **No score had been loaded.** Phase B never started.
- **No output had been written.** The occupancy matrices it computed existed only in the memory of the crashed process
  and were never written, displayed or summarised.
- **Its only output** is the per-replicate docked-frame counts that Phase A prints as it runs. These are score-free and
  say nothing about the occupancy of any individual pair.

This rule is therefore fixed before any production occupancy value has been seen, and before any correlation of any
kind has been computed from production data.

## Occupancy counts reported beside every correlation; no count threshold gates reporting — 24 October 2026

Recorded before any implementing code was written, and before any production occupancy has been written to disk,
displayed or summarised. (The first authorised attempt computed occupancy only in the memory of a process that
crashed in Phase A, before writing anything and before loading any score; see the preceding section.)

**Why.** A Spearman rho among crystal-negative pairs is defined as soon as one pair's occupancy differs from the rest.
Most crystal-negative pairs are expected to have zero occupancy, so a defined rho can rest on a handful of non-zero
pairs. Under restraint, each system's ~4,500–5,400 crystal-negative pairs take only 9–22 distinct occupancy values,
and the 6–8 A strata only 2–5 (`results/review/MD_UNDEFINED_P_REPORT.md`). On synthetic data, a single non-zero pair
among 4,970 that happens to carry the top score gives rho = 0.025 with p = 0.0005. That p is valid for the
permutation null, but it is a statement about one pair.

**What is added.** Every correlation output table — `md_occupancy_rho.csv`, `md_occupancy_strata.csv` and
`md_occupancy_docked_rho.csv` — reports, for every row, alongside `n_distinct_occupancy` and the group's n:

| Column | Counts the pairs in the row's group with |
|---|---|
| `n_nonzero_occupancy` | occupancy > 0 |
| `n_occupancy_ge_025` | occupancy >= 0.25 (the fixed meaningful-occupancy threshold) |

The counts are taken over exactly the occupancy values the row's correlation is computed on:
- **per-replicate rows:** that replicate's occupancy;
- **system rows:** the mean of the three replicates;
- **pooled rows:** the union of the systems' crystal-negative pairs;
- **strata rows:** the stratum's pairs;
- **docked rows:** the docked-frame occupancy of the contributing replicates. NaN where no replicate contributes, since
  there is then no docked occupancy to count.

**No minimum count gates reporting.** No threshold on the number of non-zero or >= 0.25 occupancies decides whether
a correlation is computed, tested or shown. Every correlation is reported with its counts. The two rules already in
force are unchanged: the undefined-correlation rule (preceding section) and the 20-pair minimum below which a distance
stratum is reported untested (decision 4).

**How such a correlation is described.** A correlation resting on very few non-zero occupancies is to be described
in the text as **uninterpretable**, not quoted as a result. It stays in the tables, with its rho, p and counts, so the
reader can see it and see why it is not interpreted.

Nothing else changes: the statistic, pair universe, window, permutation scheme, strata, thresholds, docked-frame
rules, prior control and the rows of every table are as already fixed.

## "Very few" fixed at 10: the description rule for sparse correlations — 24 October 2026

Recorded **before any implementing code was written**, and **before any production occupancy was seen**. The first
authorised attempt computed occupancy only in the memory of a process that crashed in Phase A, before writing,
displaying or summarising any of it and before loading any score (see "Undefined correlations are reported as
undefined"). No production occupancy value, count or correlation has been written, displayed or summarised since.

**The rule.** A correlation resting on **fewer than 10 non-zero occupancies** — `n_nonzero_occupancy` < 10, counted
over exactly the occupancy values the row's rho uses — is described in the text as **uninterpretable**. This fixes
the "very few" of the preceding section.

**It is a description rule only.** It is not a reporting rule and not a testing rule:
- Every such correlation is still reported in the tables, with its rho, its p (where one is computed) and all its
  counts (`n_distinct_occupancy`, `n_nonzero_occupancy`, `n_occupancy_ge_025`). **Nothing is suppressed**, and no
  permutation test is skipped because of it.
- It does not alter the existing rules. Undefined correlations still report rho and p as NaN. Distance strata under
  20 pairs are still reported untested.

**The flag.** So the rule is applied mechanically rather than by reading counts while writing, every row of
`md_occupancy_rho.csv`, `md_occupancy_strata.csv` and `md_occupancy_docked_rho.csv` carries a boolean column
`interpretable`:

> `interpretable` = rho is defined (not NaN) **and** `n_nonzero_occupancy` >= 10

It is therefore **False** for:
- a defined rho on fewer than 10 non-zero occupancies — the uninterpretable case;
- an undefined rho;
- a stratum reported untested;
- a docked row with no contributing replicate.

In the last three there is no rho to interpret, and each row's own NaN rho and counts say which case applies. The
flag changes no value in any other column. Text describing the MD results quotes only rows with
`interpretable = True` as results.

## Characterisation recomputed with periodic-boundary handling — 24 October 2026

Recorded before any fix was written or any characterisation measure recomputed.

**The defect.** `md/analysis/characterise.py` computes backbone RMSD, interface RMSD, per-loop RMSF and buried SASA
from the raw production coordinates (`prod.xtc`) without first making molecules whole or removing jumps across the
periodic box. When a chain is wrapped into a neighbouring periodic image, those coordinate-based measures see a
displaced chain and are wrong for that frame.

**How it was found.** It was found while preparing the SI figure (`results/review/FIG7_REPORT.md`, section 4). The
backbone RMSD in `md_characterisation_timeseries.csv` jumps to 15–19 A and back dozens of times in 3vxs rep1 (541
frames > 10 A, 18 separate runs), 7n1f rep3 (568 frames, 65 runs) and 7n1f rep1 (35 frames), with medians of
2.8–4.9 A between the excursions.

The decisive comparison is with `md/qc_replicate.sh`, which makes the trajectory whole (`gmx trjconv -pbc whole`,
then `-pbc nojump`) before computing backbone RMSD on frames every 10 ns:
- for **7n1f rep3 at t = 110 ns** the analysis reported **14.69 A**;
- the QC, on frames that include 110 ns, recorded a maximum of **0.45 nm** for the whole replicate (round 3).

**What it affects.** Only the descriptive characterisation measures that use coordinates directly: backbone RMSD,
interface RMSD, per-loop RMSF and buried SASA.

It does **not** affect:
- **occupancy, or any correlation computed from it** — the contact function uses minimum-image distances, verified
  by gate G2's periodic-boundary test;
- the docked-frame counts — periodic-aware;
- the interface hydrogen bonds — MDAnalysis applies the periodic box;
- any score, the primary statistic, or any table of Figure 7.

**What is done.** The four coordinate-based measures are recomputed on a trajectory made whole and jump-free by
exactly the procedure of `md/qc_replicate.sh`: `gmx convert-tpr` to the protein-only SOLU group, then
`gmx trjconv -pbc whole`, then `-pbc nojump`, on every frame rather than every nanosecond. Their definitions,
window, stride, selections and references are otherwise unchanged. The interface hydrogen bonds continue to be
computed from the original trajectory, so they are unchanged by construction; any change in them would indicate an
error. The superseded characterisation tables are kept as backups.

**This is not a second attempt at the predeclared analysis.** None of these measures sees a score, none enters the
primary statistic, and all were predeclared as descriptive. The following are **not touched**, and their hashes are
checked before and after:
- the predeclared primary statistic and every other correlation;
- the occupancy tables (`md_occupancy_pairs.csv`, `_rho`, `_strata`, `_topk`, `_check1_positive`,
  `_docked_frames`, `_docked_rho`);
- the run manifest.

The run manifest's `md_characterisation_summary_sha256_before_scores` continues to describe the superseded file.

## mdrun interruptions recorded in the QC log — 28 October 2026

Recorded after the single authorised run; it concerns only how the run history was recorded and changes no value,
table or conclusion.

The paragraph "All 12 replicates against the predeclared protocol" (24 October) lists interruptions for 5yxn rep1,
5yxn rep3, 7n1f rep2 and 3vxs rep2 only. The QC log of all twelve replicates
(`results/review/characterisation_fix/S2_qc_all12.log`, `md/qc_replicate.sh`) also records interruptions in two
replicates that paragraph does not list:

| Replicate | Queue events | mdrun sessions |
|---|---|---|
| 5yxn rep2 | 1 CRASH | 3 |
| 7n1f rep1 | 1 FAILED | 2 |

(For completeness, the log records 5yxn rep1 as 2 FAILED events over 3 sessions; the earlier text names one crash at
177 ns.)

**All 12 replicates still meet the predeclared protocol.** Both replicates above, like every other, reached the full
50,000,000 steps, wrote 10,001 frames at 20 ps with no irregular timesteps and logged 0 LINCS warnings. The protocol
allows resumption from checkpoint, which continues the same trajectory without altering its ensemble. No replicate is
excluded, re-run or re-weighted, and no occupancy, characterisation or statistical result is affected.
