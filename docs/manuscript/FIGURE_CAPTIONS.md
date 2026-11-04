# Figure captions

Final captions for the manuscript and Supporting Information figures, in figure order. **This file is the single source
for caption text**; the draft legends in the review reports point here.

- Every number below is checked against the committed tables by `scripts/58_verify_captions.py` (log:
  `results/review/captions/verify_captions.log`). The script extracts every number from every caption and fails if
  any is not traced to a table cell or to a stated method constant.
- Terminology follows the figures' own labels. Spelling is American.
- Statistics: "paired cluster bootstrap" means resampling the 69 TCR clusters (single linkage at an indel
  similarity ratio of 0.70) with replacement, pooling their complexes, and taking the mean paired difference; 95 %
  percentile confidence intervals; two-sided p floored at 1 / (number of resamples), reported as p < 0.001.
- Numbering: main Figures 1–7; Supporting Figures S1–S9 follow task order, not yet citation order.

---

## Figure 1

**Figure 1. Structural benchmark and sequence-based contact prediction with TCR-IMAP.** (A) The benchmark and the evaluation it feeds. STCRDab TCR–peptide–MHC class I complexes provide the structures,
from which the crystal contact labels and the aligned TCR and MHC sequence inputs are derived. The illustrated
complex is 1AO7: the TCR variable domains, the peptide and the MHC class I platform are shown, with TCRα in
blue, TCRβ in purple, MHC in green and the peptide in magenta; domains outside this interface view, including
β2-microglobulin and the α3 domain, are omitted. A contact is a minimum inter-residue heavy-atom distance of 4.5
Å. TCR-IMAP, the reimplemented comparators and the within-fold frequency priors are then trained and evaluated
under the five frozen random and grouped partitioning protocols, compared as paired per-complex differences and
assessed by cluster bootstrap; the primary protocol is tcr70. Protocol definitions, candidate masks, the
per-residue task and the statistical procedure are given in Methods. (B) The pairwise TCR-IMAP architecture. The
TCR side embeds 75 aligned residue rows as residue, segment and position terms and the MHC side 84 aligned
residue columns as residue and position terms; separate Transformer encoders then produce residue
representations. In the cross-attention block the TCR representations act as queries (Q) and the MHC
representations as keys (K) and values (V), and a residual connection with layer normalization updates the TCR
representations. The bilinear head combines the updated TCR representations with the MHC representations, adds a
learned pair-position term, a scale and a bias, and a sigmoid (σ) maps the result to contact probabilities
between 0 and 1. Both maps are 1AO7: the labels in (A) are its crystal contacts and the map in (B) is its stored
out-of-fold prediction under the primary tcr70 protocol, not a drawn illustration. Grey cells are residue pairs
left unresolved by the structure, which are masked in training and in scoring, and both maps keep the original
row and column order. Structural coordinates define the labels and the mask and are not input features of this
sequence-based model. The architecture shown is the pairwise variant, which predicts the complete TCR-residue ×
MHC-residue map; the predeclared primary comparison uses the separate per-residue variant, which predicts for
each TCR position whether it contacts any window residue and so matches the output format of the published
comparators (Methods). Descriptive throughout; no performance value appears in this figure.

---

## Figure 2

**Figure 2. Sequence redundancy in the 174-complex benchmark and its leakage across partitioning protocols.**
(A) V-gene usage: the eight most frequent TRAV and TRBV genes (most frequent, TRAV12 in 61 complexes and TRBV6 in
49). (B) Cumulative distributions of the TCR similarity ratio over all 15,051 complex pairs and for each
complex's nearest neighbor. A ratio ≥ 0.90 holds for 1.5 % of pairs but for the nearest neighbor of 67.8 % of
complexes, and 100 of 174 complexes (57.5 %) have an exact TCR twin. (C) Complexes by cluster size at 70, 80 and 90 %
similarity (69, 82 and 89 clusters). (D) For each test complex, the maximum similarity ratio to any complex
in its training folds, under each protocol. Under random partitioning the median test complex has a training complex at
0.983 and 60.3 % have one at ≥ 0.90; under tcr70 the maximum is 0.696, and no test complex reaches 0.90. (E) Duplicate
families (identical CDR3α–CDR3β pair and peptide) split between training and test, mean per fold: random 7.4, intersection (legacy)
0.67, allele 1.25, temporal 0 and tcr70 0. Complexes evaluated under each protocol are given below each bar (random
174, intersection (legacy) 30, allele 64, temporal 60 and tcr70 174 of 174). The 108 complexes without an allele
annotation fall outside every allele and intersection test fold, and one allele fold of two complexes is below the
minimum fold size. The partition labeled intersection is a legacy intersection-selected partition, not a double
holdout: its held-out groups remain in training (Methods). The
similarity ratio is computed on the concatenated CDR1, CDR2 and CDR3 slot sequences of both chains, one residue per
IMGT position (IMGT 27–38, 56–65 and 105–124, the last including J-region framework positions 118–124), the sequences
used to cluster the complexes and build the folds. 6UZ1 is excluded from the benchmark on structural-scope grounds.
Descriptive; no test is performed. n = 174 complexes.

---

## Figure 3

**Figure 3. Apparent performance depends on train/test redundancy.** (A) Pairwise macro-AUPRC of TCR-IMAP and the
pair-frequency prior when every training complex whose TCR similarity ratio to any test complex of its
fold reaches the threshold is removed (none, no restriction). One fixed random five-fold partition is used throughout,
and all 174 complexes are evaluated at every threshold. From no restriction to 0.70, TCR-IMAP falls by 0.3042 (0.6313
to 0.3271) and the prior by 0.0566 (0.3563 to 0.2997). (B) TCR-IMAP − prior with 95 % confidence intervals. The
advantage decreases at every step, from +0.2750 (p < 0.001) without restriction to +0.0274 at 0.70 ([+0.0011,
+0.0542], p = 0.041). Filled markers, CI excludes zero. Intervals and p values use 20,000 bootstrap resamples at
seed 42; the 2,000-resample values first recorded for this stage gave [−0.0001, +0.0543], p = 0.052 at this threshold
from the same point estimate, a Monte Carlo difference in the bootstrap, and are retained in the Supporting
Information. (C) Mean absolute error (MAE) of the
docking angles predicted by AngleNet, a separate sequence-based regression network trained from scratch in each fold
(not TCR-IMAP), under random, allele-holdout and TRBV-holdout partitions: three training seeds each (points) and their
mean (line). MAE 5.40–6.08° (random), 11.18–12.01° (allele-annotation holdout) and 12.50–13.43° (TRBV holdout); dashed line, the
median-predictor baseline (11.0°). Pearson r across seeds is printed for each partition: 0.71 to 0.78, 0.13 to 0.20 and
0.11 to 0.20. (D) Docking-angle R² for the same runs: +0.42 to +0.56 under random partitioning, −0.38 to −0.22 under
allele-annotation holdout and −0.47 to −0.39 under TRBV holdout; dashed line, R² = 0 (the mean predictor). Statistics: (B) paired
cluster bootstrap over 69 TCR clusters, 20,000 resamples, seed 42, recomputed from per-complex arrays recovered under
a validation gate that reproduced the frozen aggregates exactly; the superseded 2,000-resample values are retained in
the Supporting Information. (E) Training-size-matched control at the
three thresholds for which it was run. Squares, the similarity-filtered advantage; small points, the five fixed
size-matched random draws; large circles, the principal size-matched summary, the mean over draws of the per-complex
model-minus-prior difference. Restricting similarity also shrinks training, from 139–140 complexes per fold without
restriction to 91–96, 86–93 and 69–78 at thresholds of 0.95, 0.85 and 0.70; each control draws that same number of
training complexes at random from the unrestricted pool, with no similarity filter, on the same test complexes. The
size-matched advantage is +0.2400 [+0.1827, +0.2832], +0.2364 [+0.1748, +0.2822] and +0.2296 [+0.1687, +0.2747]
against the filtered +0.0777, +0.0428 and +0.0274; the paired contrast between the two regimes, computed per complex
on the same test complexes under the same resamples, is −0.1623 [−0.2090, −0.1163], −0.1936 [−0.2434, −0.1364] and
−0.2021 [−0.2455, −0.1471], all p < 0.001 (printed beside each pair). Dotted line, the unrestricted advantage
(+0.2750). Filled markers, 95 % CI excludes zero. Error bars are 95 % cluster-bootstrap intervals (20,000 resamples,
seed 42); the five draws are plotted as points because their spread is a range over draws, not a confidence interval.
Statistics: (B) and (E) paired cluster bootstrap over 69 TCR clusters, 20,000 resamples, seed 42. n = 174 complexes (A, B, E); 172 of 174 (C, D; 3vxu and 3w0w have no docking angle).

---

## Figure 4

**Figure 4. TCR-IMAP against published comparators.** (A) Per-residue macro-AUPRC of TCR-IMAP, Struct-TCS, Seq-TCS,
the within-fold position prior and nearest-neighbor (NN) transfer under five partitioning protocols, ordered by
partitioning protocols, with the number of complexes evaluated (random 174, intersection (legacy) 30, allele 64,
temporal 60, tcr70 174). The protocols are not ordered by strictness: they probe different distribution shifts on
different subsets. Each protocol is annotated with the TCR-IMAP and Struct-TCS means: random 0.806 and 0.809,
intersection (legacy) 0.800 and 0.730, allele
0.776 and 0.726, temporal 0.694 and 0.691, tcr70 0.728 and 0.697. (B) Paired difference TCR-IMAP − Struct-TCS per
protocol with 95 % confidence intervals: random −0.0030 (p = 0.798), intersection (legacy) +0.0696 (p = 0.015), allele +0.0507
(p = 0.001), temporal +0.0030 (p = 0.874) and tcr70, the predeclared primary comparison, +0.0312 [−0.0049, +0.0650]
(p = 0.088), which is not significant. (C, D) Pairwise contact prediction under tcr70, one point per complex; dashed,
identity line. (C) TCR-IMAP against the within-fold pair-frequency prior: TCR-IMAP higher on 112 of 174 complexes
(64.4 %), pooled Δ = +0.0507 [+0.0215, +0.0841], p < 0.001. (D) TCR-IMAP against TEIM-Res: higher on 143 of 174
(82.2 %), pooled Δ = +0.0873 [+0.0592, +0.1149], p < 0.001. (E) Cluster-macro sensitivity, tcr70: the three primary
comparisons with complexes weighted equally (complex-macro, circles, n = 174; the predeclared weighting) and with the
69 TCR clusters weighted equally (cluster-macro, squares). Per-residue vs Struct-TCS, p = 0.088 and p = 0.063;
per-residue vs the position prior, p = 0.546 and p = 0.174; pairwise vs the pair-frequency prior, p < 0.001 and
p = 0.001. The two weightings agree in every comparison. (F) The 88 structures shared with the comparator's dataset,
at its 4.5 Å and 6 Å cutoffs, under a single random split of the comparator's design (n = 22 test complexes) and under
tcr70 (n = 88). Only the structure subset and cutoff are shared; candidate positions, TCR representation, task pooling
and validation differ, so this is not a reproduction of the comparator's own evaluation. Under tcr70, TCR-IMAP leads
at both cutoffs (+0.1057 [+0.0705, +0.1357] and +0.0964 [+0.0684, +0.1235], both p < 0.001). Under the single random
split the difference is not significant at 4.5 Å (+0.0269 [−0.0102, +0.0707], p = 0.171) and reverses at 6 Å
(−0.0949 [−0.1335, −0.0517], p < 0.001). All per-residue models in A, B, E and F are selected on cluster-disjoint
validation. Statistics: paired cluster bootstrap over the 69 TCR clusters, 20,000 resamples, seed 42, in every panel;
filled markers, 95 % CI excludes zero.

---

## Figure 5

**Figure 5. Observed structural contact footprint on the MHC α1/α2 platform.** (A) Contact frequency (contacts /
complexes in which both residues are resolved; minimum heavy-atom distance ≤ 4.5 Å) for 75 TCR rows (68 occupied IMGT
positions and 7 CDR3 insertion positions) × the 84-column α1/α2 window. The window is IMGT 50–90 and 1050–1090 plus the
insertion residues 1061A and 1072A (41 α1 and 43 α2 columns); α2 positions are shown without the 1000 offset. The strip
marks DSSP helix (dark) and non-helical boundary (light) columns. Boxes mark the helix each loop prefers (see E); solid
outline, germline loops; dashed, CDR3; color, TCR chain. Three CDR3β insertion rows resolved in fewer than five
complexes (111B and 112C, one complex each; 112B, three) are hatched and kept off the color scale, because a frequency
over one to three complexes is not a population estimate. With these rows excluded, the most frequent pair is CDR1β 36
× α1 76 (0.667). Square-root color scale. (B) The MHC platform of 1AO7 colored by each position's contact frequency
summed over the TCR rows resolved in at least five complexes, relative to its own maximum (not panel A's scale).
Sticks mark the 14 positions above 0.30 of the maximum: 62, 65, 66, 68, 69, 72, 73, 76, 1061A, 1062, 1065, 1066, 1069
and 1073. Residues are selected by IMGT number and insertion code. (C) The 1AO7 TCR–pMHC complex with its CDR loops.
(D) Contact rate (%) by loop and MHC region over the 182-position G-domain. The β-sheet floor carries 0 contacts in
595,952 resolved pairs; the 84-column window captures 5,281 of the 5,294 G-domain contacts (99.75 %). (E) Per-cluster
α1 − α2 contact-rate difference with 95 % bootstrap intervals over the 69 TCR clusters. α-chain germline loops prefer
α2 and β-chain germline loops α1; CDR3α prefers α1 and CDR3β α2. CDR2β contacts no α2-helix position (0 of 40,560
resolved pairs). Statistics: (E) cluster bootstrap, 2,000 resamples, seed 42; the p values shown are unadjusted
(p < 0.001 for five loops, p = 0.002 for CDR3β), and after Holm correction across the six loops every loop has
p = 0.003; filled markers, CI excludes zero. Descriptive: no model and no data split. α3 is excluded. n = 174
complexes.

---

## Figure 6

**Figure 6. Sequence dependence of TCR-IMAP, and the complexes on which it gains and loses most against the
pair-frequency prior (tcr70, pairwise task).** (A) Change in pairwise macro-AUPRC when one input or architectural
component is removed (full − ablated; full TCR-IMAP 0.3598). Bars are Bonferroni-simultaneous 99.375 % intervals
across the eight ablations; raw and Holm-adjusted p values are given beside each row. Filled markers, simultaneous
interval excludes zero, which coincides with Holm-adjusted p < 0.05 in every row. Three ablations lowered performance
after correction: TCR residues shuffled within loops (+0.0503, Holm p < 0.001), CDR1/2 removed (+0.0461, Holm
p < 0.001) and MHC residues shuffled across positions (+0.0283, Holm p = 0.014). Masking TCR residue identity (+0.0337,
Holm p = 0.104) or MHC residue identity (+0.0180, Holm p = 0.160) did not survive correction. Removing cross-attention
(−0.0008) or replacing the bilinear head with a dot product (−0.0039) had no detectable effect, so neither component is
shown to be necessary. (B) Per-complex difference between TCR-IMAP and the within-fold pair-frequency prior, ranked:
TCR-IMAP is higher on 112 of 174 complexes (64.4 %); pooled Δ = +0.0507 [+0.0215, +0.0841], p < 0.001 (dashed line).
(C, D) Observed contacts (left) and TCR-IMAP's k highest-scoring residue pairs (right), where k is the complex's number
of observed contacts, from out-of-fold tcr70 predictions. View onto the MHC α1/α2 platform from the TCR, α1 helix at
the top, peptide N-terminus at the left. MHC surface residues are colored by the TCR chain contacting them (sky blue,
TCRα; vermillion, TCRβ); CDR loops are drawn as tubes, with contact residues as spheres at Cα. MHC residues contacted
by both chains are drawn in both colors. (C) 3QFJ, the largest gain over the prior (pairwise AUPRC 0.566 vs 0.287,
Δ = +0.2793; k = 30; 16 of 30 predicted pairs observed; MHC sites shared, 13 of 17). (D) 3QEQ, the largest loss to it
(0.293 vs 0.452, Δ = −0.1583; k = 33; 10 of 33 observed; 9 of 14). Overlap counts equal precision at k for these two
complexes and are illustrative, not a performance estimate. Contacts are residue pairs within 4.5 Å (minimum
heavy-atom distance) on the 84-column α1/α2 window. Statistics: paired cluster bootstrap over 69 TCR clusters, 20,000
resamples, seed 42. n = 174 complexes.

---

## Figure 7

**Figure 7. Solution-phase contact of crystal-negative residue pairs, against TCR-IMAP scores and the pair-frequency
prior.** Four TCR–pMHC complexes, selected before any simulation for a high (5YXN, 3VXS) or low (7N1F, 8EO8)
per-complex gain of TCR-IMAP over the prior, were simulated in three independent 200 ns replicates each. The occupancy
of a residue pair is the fraction of the 9,001 frames between 20 and 200 ns in which its minimum heavy-atom distance
is ≤ 4.5 Å, averaged over the three replicates. Scores are TCR-IMAP's out-of-fold tcr70 pairwise predictions and the
within-fold pair-frequency prior. (A) Spearman ρ between score and occupancy over crystal-negative resolved pairs,
grouped by the pair's minimum heavy-atom distance in the crystal structure (crystal-frame distance stratum), with the
four systems pooled. Bars, pooled ρ; dots, single systems, filled when the approximate within-system permutation
p < 0.05 (10,000 permutations of the scores among a system's crystal-negative pairs, two-sided; approximate because
pairs within a system are not independent). Beneath each group, its pair count and number of non-zero occupancies. Over
all pairs (n = 19,957, 381 non-zero), occupancy correlates more strongly with the prior than with TCR-IMAP (ρ = 0.44 vs
0.22), in every system and every replicate. The difference arises in the > 8 Å stratum (0.21 vs 0.10), which holds
19,414 of the 19,957 pairs, of which 94 have any contact and none reaches 0.25. The two scores are similar in the 6–8 Å
stratum (0.31 vs 0.30; 423 pairs). In the 4.5–6 Å stratum (120 pairs) TCR-IMAP is numerically higher (0.24 vs 0.18),
but no test of the difference between scores was made, so no stratum shows either score to rank transient contacts
better. No correlation rests on fewer than 10 non-zero occupancies. Across these pairs the two scores are only
moderately correlated with each other (Spearman ρ = 0.33–0.39 per system, 0.36 pooled). (B) Score against occupancy for
every crystal-negative pair, colored by stratum. Non-zero occupancies on a log scale; the 19,576 pairs with zero
occupancy in the strip below. Dashed line, occupancy 0.25. (C) Fraction of the k highest-scoring crystal-negative pairs
with occupancy ≥ 0.25 (colored segments, labeled with counts) against the fraction among all crystal-negative pairs
(base rate, dotted line), per system. k is the system's crystal-positive count (22–36). Both scores enrich
high-occupancy pairs 28- to 88-fold over the base rate, on 1–7 pairs per system; † marks a tie in score at rank k.
(D) Occupancy of the crystal-positive pairs (check 1) per replicate and per system (all): line, mean; ×, median. Above
each replicate, the percentage of frames that keep at least 50 % of the starting TCR–pMHC contacts, read two ways:
docked % (count), the frame has at least half as many TCR–pMHC heavy-atom contacts as at the start; docked %
(identity), at least half of the specific starting contact pairs are present. The TCR of 5YXN disengages in two of its
three replicates (5 % and 3 % of frames docked by count), and its crystal contacts are the least occupied (system mean
0.33). (E) Backbone RMSD (protein Cα, from the energy-minimized structure) over 20–200 ns for each replicate; dashed
line, system mean (5YXN 4.25, 7N1F 3.98, 3VXS 3.12 and 8EO8 3.80 Å); the docked % (count) of each replicate is printed
above. Backbone RMSD stays between 1.5 and 7.7 Å in every replicate, including the 5YXN replicates in which the TCR
disengages: that is a change of interface, not a loss of fold (Figure S3). Backbone RMSD is shown as evidence of
simulation stability; it is descriptive and is not used to relate stability to prediction quality.

---

## Figure S1

**Figure S1. Effect of adding the peptide sequence to TCR-IMAP: no detectable effect.** Difference in pairwise
macro-AUPRC under tcr70 between TCR-IMAP trained with and without the peptide sequence, with 95 % confidence intervals,
in four scopes: all positions (+0.0100 [−0.0058, +0.0258], p = 0.196); all positions with the peptide sequence shuffled,
a control (+0.0071 [−0.0072, +0.0222], p = 0.308); CDR1/2 rows only (+0.0107 [−0.0038, +0.0267], p = 0.133); and CDR3
rows only (+0.0082 [−0.0137, +0.0293], p = 0.416). Every 95 % CI includes zero (hollow markers). The macro-AUPRC without
→ with peptide is given beside each row. Statistics: paired cluster bootstrap over 69 TCR clusters, 2,000 resamples,
recorded values (no per-complex arrays are stored for this analysis). n = 174 complexes.

---

## Figure S2

**Figure S2. Per-complex per-residue AUPRC under each partitioning protocol.** Distribution over test complexes of the
per-residue AUPRC of TCR-IMAP, Struct-TCS, Seq-TCS, the within-fold position prior and nearest-neighbor (NN) transfer,
under random (n = 174), intersection/legacy (n = 30), allele (n = 64), temporal (n = 60) and tcr70 (n = 174)
partitioning. Boxes
show the median and interquartile range (IQR), whiskers extend to 1.5 × IQR, and points are individual complexes.
The means are those of Figure 4A. Per-complex AUPRC tracks a complex's number of contacts and is not a measure of
difficulty; differences between methods are assessed with paired tests (Figure 4B). Descriptive; no test.

---

## Figure S3

**Figure S3. Interface characterization of the four simulated TCR–pMHC complexes.** Columns, systems (5YXN, 7N1F,
3VXS, 8EO8); colors, the three independent replicates. Every measure covers the analyzed window, 20–200 ns of each
200 ns production run, on trajectories made whole and free of periodic jumps. Backbone RMSD is shown in Figure 7E.
(A) Interface RMSD: RMSD of the Cα atoms of the canonical CDR loops (IMGT 27–38, 56–65 and 105–117) after a
least-squares fit on the MHC α1/α2 Cα atoms, relative to the crystal pose; every frame, 20 ps apart (9,001 per
replicate). (B) Buried solvent-accessible surface area (SASA) of the TCR–pMHC interface, SASA(TCR) + SASA(pMHC) −
SASA(complex), Shrake–Rupley with a 1.4 Å probe, sampled every 1 ns (181 points per replicate), which is why these
traces are smoother than those in A. (C) Per-loop Cα RMSF after the same fit over all 9,001 frames, averaged within
each CDR loop (loop labels colored by chain). (D) Interface hydrogen bonds per frame (donor–acceptor ≤ 3.5 Å,
D–H···A ≥ 150°, every frame), between the TCR and MHC α1/α2 and between the TCR and the peptide: bars, mean over the
window; whiskers, frame-to-frame SD (lower whisker clipped at zero). 5YXN's disengagement appears jointly in these
measures: in replicate 1 the interface RMSD exceeds 8 Å in more than a quarter of the frames (maximum 10.7 Å;
mean 6.85 Å, the highest of all twelve replicates);
buried SASA falls in all three replicates, most in replicate 2 (from about 2,010 to about 945 Å² between the first
and the last 10 ns), and 5YXN has the lowest mean of the four systems (1,316 Å²); CDR3β RMSF is elevated in
replicates 1 and 2 (3.18 and 2.02 Å), where the β chain lost its contacts, and not in replicate 3 (0.94 Å), where the
α chain disengaged instead. These measures are descriptive: with four systems they characterize these simulations
only, and they are not used to relate flexibility, burial or hydrogen bonding to prediction quality.

---

## Figure S4

**Figure S4. Sensitivity of the per-residue benchmark to the validation split used for model selection.**
(A) Per-residue macro-AUPRC of TCR-IMAP, Struct-TCS and Seq-TCS under each partitioning protocol, with early stopping
on a random 15 % of the training complexes (random split validation, triangles; superseded) and on whole TCR clusters
(cluster-disjoint validation, circles; reported). The value above each protocol is the change in TCR-IMAP's
macro-AUPRC from the random to the cluster-disjoint split: random −0.0541, intersection (legacy) −0.0252, allele −0.0484, temporal
−0.0122 and tcr70 −0.0165. The position prior and NN transfer use no validation subset and are identical under both
splits (not shown). (B) Under tcr70, TCR-IMAP − Struct-TCS, − Seq-TCS and − the position prior with 95 % confidence
intervals under each split. The predeclared primary comparison, against Struct-TCS, is +0.0504 [+0.0188, +0.0828]
(p = 0.001) under the random split and +0.0312 [−0.0049, +0.0650] (p = 0.088) under the cluster-disjoint split; against
Seq-TCS, +0.0897 and +0.0799 (both p < 0.001); against the position prior, +0.0250 (p = 0.060) and +0.0085
(p = 0.546). (C) Number of TCR clusters with complexes on both sides of model selection, in each of the 18 folds under
the random split: 5–15 (median 10.5). Under the cluster-disjoint split the number is 0 in every fold, by construction.
Statistics: (B) paired cluster bootstrap over 69 TCR clusters, 20,000 resamples, seed 42, recomputed from the stored
per-complex AUPRC arrays; filled markers, 95 % CI excludes zero. n = 174 complexes (random, tcr70), 30 (intersection),
64 (allele), 60 (temporal).

---

## Figure S5

**Figure S5. Sensitivity of the three primary tcr70 comparisons to cluster weighting.** Difference in macro-AUPRC
between TCR-IMAP and the comparator with every test complex weighted equally (complex-macro, circles; the predeclared
primary weighting) and with each of the 69 TCR clusters weighted equally, complexes being averaged within their
cluster first (cluster-macro, squares). Per-residue vs Struct-TCS: +0.0312 [−0.0049, +0.0650], p = 0.088 (complex-macro)
and +0.0258 [−0.0016, +0.0518], p = 0.063 (cluster-macro). Per-residue vs the position prior: +0.0085 [−0.0162,
+0.0371], p = 0.546 and +0.0164 [−0.0070, +0.0400], p = 0.174. Pairwise vs the pair-frequency prior: +0.0507 [+0.0215,
+0.0841], p < 0.001 and +0.0360 [+0.0148, +0.0579], p = 0.001. Complex-macro is the primary analysis, and the two
weightings agree in every comparison: neither per-residue comparison is significant under either weighting, and the
pairwise comparison is significant under both. Right, the macro-AUPRC of TCR-IMAP and the comparator under each
weighting. Filled markers, 95 % CI excludes zero. Statistics: paired cluster bootstrap over 69 TCR clusters, 20,000
resamples, seed 42, recomputed from the stored per-complex arrays. n = 174 complexes, 69 clusters.

---

## Figure S6

**Figure S6. The comparator subset, per training seed.** The 88 structures shared with the comparator's dataset, at its
4.5 Å and 6 Å contact cutoffs, under a single random split of the kind the comparator used (n = 22 test complexes) and
under tcr70 (n = 88). (A) TCR-IMAP − Struct-TCS per-residue macro-AUPRC. Circles, the reported three-seed ensemble
(predictions averaged, then scored) with model selection on cluster-disjoint validation, with 95 % CI and p; crosses,
the three training seeds individually; triangles, the same ensemble with model selection on a random validation subset
(superseded). (B) Macro-AUPRC of each method: ensemble (circles) and single seeds (crosses). Under cluster-disjoint
validation the 6 Å single-random-split cell reverses sign, from +0.0265 [+0.0008, +0.0506] to −0.0949 [−0.1335,
−0.0517] (p < 0.001), and all three seeds agree on the sign (−0.060, −0.120 and −0.087); the seeds agree in sign in
every cell. Under tcr70, TCR-IMAP leads at both cutoffs (+0.1057 and +0.0964, both p < 0.001). The single-random-split
rows are reported for comparability with the published comparator setting, not as an evaluation we endorse: each rests
on one fold of 22 test complexes, and only the structure subset and cutoff are shared with the comparator's own
evaluation. Filled markers, 95 % CI excludes zero. Statistics: paired cluster bootstrap over TCR clusters, 20,000
resamples, seed 42; the reported ensemble was recomputed from the stored per-complex arrays.

---

## Figure S7

**Figure S7. Legacy results on 28 MHC positions, shown to document the correction.** An earlier version of this
analysis used 28 MHC positions that were selected from contact frequency pooled over all 174 complexes before the data
were split into training and test folds. That is outcome-informed selection: the candidate set was chosen with
knowledge of the test complexes' contacts. It was replaced by the structurally defined 84-column α1/α2 window used
throughout this paper (p84v2), which does not depend on observed contacts. The legacy results (diamonds) are superseded
and are shown only to document what the correction changed; the p84v2 values (circles) are those reported in the main
text. Five headline comparisons, TCR-IMAP − comparator in macro-AUPRC, with 95 % CI and p: pairwise tcr70 vs the
pair-frequency prior, +0.0398 (p = 0.010) → +0.0507 (p < 0.001); per-residue tcr70 vs Struct-TCS, +0.0273 (p = 0.001)
→ +0.0312 (p = 0.088); per-residue tcr70 vs the position prior, +0.0256 (p = 0.048) → +0.0085 (p = 0.546); identity
sweep at < 0.70 vs the pair-frequency prior, +0.0127 (p = 0.429) → +0.0274 (p = 0.052); peptide input, CDR3 rows,
+0.0380 (p < 0.001) → +0.0082 (p = 0.416). The per-residue legacy models were also selected on a random rather than a
cluster-disjoint validation subset, so their change combines both corrections. Legacy values are as originally
recorded (2,000 resamples); p84v2 values use 20,000 resamples (seed 42) where per-complex arrays exist, and the recorded
2,000-resample values for the identity sweep and the peptide comparison. Filled markers, 95 % CI excludes zero.
n = 174 complexes.

---

## Figure S8

**Figure S8. Per-complex pairwise AUPRC under tcr70, grouped by MHC allele and by TRBV.** Points, individual
complexes; lines, group means; TCR-IMAP (blue) and the within-fold pair-frequency prior (gray); n above each group.
(A) By MHC allele. Allele annotation was parsed from the structure records and is missing for 108 of 174 complexes,
shown as one separate group. (B) By TRBV subgroup. Groups with fewer than three complexes are omitted (A: 11 alleles,
13 complexes; B: 6 groups, 8 complexes, including one complex without a TRBV annotation). This figure is descriptive:
it is not a subgroup analysis and no test is performed. Groups are small and not independent of the TCR clusters used
for evaluation, and per-complex AUPRC tracks a complex's number of contacts, so it is not a measure of difficulty.
n = 174 complexes.

---

## Figure S9

**Figure S9. Composition of the identity-threshold sweep.** The sweep of Figure 3A,B removes, from each training fold,
every complex whose TCR similarity ratio to any test complex of that fold reaches the threshold, on a fixed
five-fold partition. (A) TCR-IMAP − pair-frequency prior, pairwise macro-AUPRC, with 95 % CI (as Figure 3B; recorded
values, paired cluster bootstrap, 2,000 resamples); filled markers, CI excludes zero; at 0.70, +0.0274, p = 0.052.
(B) Complexes per fold. The test composition does not change: all 174 complexes are evaluated at every threshold, 34–35
per fold. The training composition does change, which is what the sweep manipulates: from 139–140 complexes per fold
without restriction to 69–78 at 0.70 (crosses, individual folds; line, mean). (C) The train/test similarity achieved:
for each test complex, its maximum similarity to the training complexes its fold kept, summarized as the maximum
(triangles) and mean (circles) over the 174 test complexes; dotted line, the threshold. At 0.70 the achieved maximum is
0.699 and the mean 0.649. The dash-dotted line is the mean_max_sim value recorded with the sweep, which is not the
achieved similarity: it takes each test complex's maximum similarity to any other complex below the threshold,
including complexes in the same test fold and training complexes the filter removed, and is therefore an upper bound on
the achieved mean. n = 174 complexes.
