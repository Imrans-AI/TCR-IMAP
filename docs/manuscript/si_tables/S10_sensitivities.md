# Table S10. Common-test-set and target-definition sensitivities

**Table S10.** Two sensitivities computed entirely from stored prediction and per-complex arrays. Nothing was
retrained for either.

## S10a. The comparator subset on a common test set

The reported comparator-subset contrast sets a single random split of 22 test complexes against tcr70 with 88, so it
confounds the partitioning protocol with the set of complexes scored. The 22 are a subset of the 88, so the tcr70
evaluation can be restricted to exactly those 22. The 88-complex rows are preserved and remain the reported result.
Paired cluster bootstrap over clu70, 20,000 resamples, seed 42.

| cutoff | comparator | setting | n | clusters | TCR-IMAP | comparator AUPRC | Δ | lo | hi | p |
|---|---|---|---|---|---|---|---|---|---|---|
| 4.5 A | Struct-TCS | single random split, n = 22 | 22 | 17 | 0.8108 | 0.7839 | 0.0269 | −0.0102 | 0.0707 | 0.1705 |
| 4.5 A | Struct-TCS | tcr70, full subset, n = 88 | 88 | 40 | 0.7297 | 0.6241 | 0.1057 | 0.0705 | 0.1357 | 0.0001 |
| 4.5 A | Struct-TCS | tcr70, SAME 22 complexes | 22 | 17 | 0.7361 | 0.6449 | 0.0912 | 0.0298 | 0.1514 | 0.0011 |
| 4.5 A | Seq-TCS | single random split, n = 22 | 22 | 17 | 0.8108 | 0.7375 | 0.0733 | 0.0308 | 0.1155 | 0.0003 |
| 4.5 A | Seq-TCS | tcr70, full subset, n = 88 | 88 | 40 | 0.7297 | 0.6023 | 0.1274 | 0.0884 | 0.1618 | 0.0001 |
| 4.5 A | Seq-TCS | tcr70, SAME 22 complexes | 22 | 17 | 0.7361 | 0.6093 | 0.1267 | 0.0658 | 0.1855 | 0.0001 |
| 6 A | Struct-TCS | single random split, n = 22 | 22 | 17 | 0.7621 | 0.8570 | −0.0949 | −0.1335 | −0.0517 | 0.0001 |
| 6 A | Struct-TCS | tcr70, full subset, n = 88 | 88 | 40 | 0.8166 | 0.7202 | 0.0964 | 0.0684 | 0.1235 | 0.0001 |
| 6 A | Struct-TCS | tcr70, SAME 22 complexes | 22 | 17 | 0.8220 | 0.7485 | 0.0735 | 0.0172 | 0.1233 | 0.0098 |
| 6 A | Seq-TCS | single random split, n = 22 | 22 | 17 | 0.7621 | 0.8333 | −0.0712 | −0.1195 | −0.0215 | 0.0036 |
| 6 A | Seq-TCS | tcr70, full subset, n = 88 | 88 | 40 | 0.8166 | 0.7392 | 0.0774 | 0.0545 | 0.1045 | 0.0001 |
| 6 A | Seq-TCS | tcr70, SAME 22 complexes | 22 | 17 | 0.8220 | 0.7611 | 0.0609 | 0.0152 | 0.1153 | 0.0067 |

At 6 Å the sign of the difference from Struct-TCS **reverses between protocols on identical test complexes** (−0.0949 under the single random split against +0.0735 under tcr70).
That locates the instability in the partitioning design rather than in the choice of test complexes.
Seed-to-seed variation within a fixed split (Figure S6) does not measure variation between splits and is not used to infer it.

## S10b. Canonical-CDR-only scoring, FR4 excluded

The 75-row TCR representation keeps 14 J-region framework rows (IMGT 118–124 of both chains).
Those rows carry **0 of the 5,281 contacts** at 4.5 Å but **204,624 of the 893,928 resolved pairs**, 22.9 % of the evaluated negatives and none of the positives.
Rescoring the stored predictions over canonical CDR rows only (IMGT ≤ 117) tests how much each comparison depends on them.
Paired cluster bootstrap over clu70, 20,000 resamples, seed 42. The *TCR-IMAP* and *comparator* columns are macro-AUPRC.

| evaluation | comparison | scope | n | TCR-IMAP | comparator | Δ | lo | hi | p |
|---|---|---|---|---|---|---|---|---|---|
| pairwise, random | TCR-IMAP vs pair-frequency prior | all rows (reported) | 174 | 0.6230 | 0.3559 | 0.2671 | 0.2069 | 0.3126 | 0.0001 |
| pairwise, random | TCR-IMAP vs pair-frequency prior | canonical CDRs only | 174 | 0.6230 | 0.3560 | 0.2671 | 0.2069 | 0.3126 | 0.0001 |
| pairwise, allele | TCR-IMAP vs pair-frequency prior | all rows (reported) | 64 | 0.4721 | 0.3306 | 0.1415 | 0.0577 | 0.2241 | 0.0001 |
| pairwise, allele | TCR-IMAP vs pair-frequency prior | canonical CDRs only | 64 | 0.4722 | 0.3307 | 0.1415 | 0.0577 | 0.2242 | 0.0001 |
| pairwise, intersection (legacy) | TCR-IMAP vs pair-frequency prior | all rows (reported) | 30 | 0.4406 | 0.3182 | 0.1223 | 0.0508 | 0.2021 | 0.0002 |
| pairwise, intersection (legacy) | TCR-IMAP vs pair-frequency prior | canonical CDRs only | 30 | 0.4406 | 0.3183 | 0.1223 | 0.0507 | 0.2021 | 0.0002 |
| pairwise, temporal | TCR-IMAP vs pair-frequency prior | all rows (reported) | 60 | 0.3200 | 0.2892 | 0.0308 | 0.0012 | 0.0663 | 0.0414 |
| pairwise, temporal | TCR-IMAP vs pair-frequency prior | canonical CDRs only | 60 | 0.3203 | 0.2893 | 0.0310 | 0.0014 | 0.0665 | 0.0400 |
| pairwise, tcr70 | TCR-IMAP vs pair-frequency prior | all rows (reported) | 174 | 0.3598 | 0.3090 | 0.0507 | 0.0215 | 0.0841 | 0.0001 |
| pairwise, tcr70 | TCR-IMAP vs pair-frequency prior | canonical CDRs only | 174 | 0.3599 | 0.3091 | 0.0508 | 0.0215 | 0.0841 | 0.0001 |
| pairwise, tcr70 | TCR-IMAP vs TEIM-Res | all rows (reported) | 174 | 0.3598 | 0.2725 | 0.0873 | 0.0592 | 0.1149 | 0.0001 |
| pairwise, tcr70 | TCR-IMAP vs TEIM-Res | canonical CDRs only | 174 | 0.3599 | 0.2740 | 0.0859 | 0.0580 | 0.1133 | 0.0001 |
| per-residue, tcr70, 88-complex comparator subset | TCR-IMAP vs Struct-TCS | all rows (reported) | 88 | 0.7297 | 0.6241 | 0.1057 | 0.0705 | 0.1357 | 0.0001 |
| per-residue, tcr70, 88-complex comparator subset | TCR-IMAP vs Struct-TCS | canonical CDRs only | 88 | 0.7307 | 0.6248 | 0.1059 | 0.0706 | 0.1357 | 0.0001 |
| per-residue, tcr70, 88-complex comparator subset | TCR-IMAP vs Seq-TCS | all rows (reported) | 88 | 0.7297 | 0.6023 | 0.1274 | 0.0884 | 0.1618 | 0.0001 |
| per-residue, tcr70, 88-complex comparator subset | TCR-IMAP vs Seq-TCS | canonical CDRs only | 88 | 0.7307 | 0.6038 | 0.1269 | 0.0884 | 0.1609 | 0.0001 |
| sweep 0.95, size-matched random | TCR-IMAP vs pair-frequency prior | all rows (reported) | 174 | 0.6087 | 0.3570 | 0.2517 | 0.1906 | 0.2979 | 0.0001 |
| sweep 0.95, size-matched random | TCR-IMAP vs pair-frequency prior | canonical CDRs only | 174 | 0.6087 | 0.3570 | 0.2517 | 0.1906 | 0.2979 | 0.0001 |
| sweep 0.95, similarity-filtered | TCR-IMAP vs pair-frequency prior | all rows (reported) | 174 | 0.3949 | 0.3172 | 0.0777 | 0.0336 | 0.1166 | 0.0004 |
| sweep 0.95, similarity-filtered | TCR-IMAP vs pair-frequency prior | canonical CDRs only | 174 | 0.3950 | 0.3173 | 0.0777 | 0.0337 | 0.1167 | 0.0004 |
| sweep 0.85, size-matched random | TCR-IMAP vs pair-frequency prior | all rows (reported) | 174 | 0.6099 | 0.3575 | 0.2524 | 0.1825 | 0.3038 | 0.0001 |
| sweep 0.85, size-matched random | TCR-IMAP vs pair-frequency prior | canonical CDRs only | 174 | 0.6100 | 0.3576 | 0.2525 | 0.1826 | 0.3039 | 0.0001 |
| sweep 0.85, similarity-filtered | TCR-IMAP vs pair-frequency prior | all rows (reported) | 174 | 0.3526 | 0.3098 | 0.0428 | 0.0101 | 0.0747 | 0.0087 |
| sweep 0.85, similarity-filtered | TCR-IMAP vs pair-frequency prior | canonical CDRs only | 174 | 0.3528 | 0.3099 | 0.0428 | 0.0101 | 0.0747 | 0.0087 |
| sweep 0.7, size-matched random | TCR-IMAP vs pair-frequency prior | all rows (reported) | 174 | 0.5995 | 0.3565 | 0.2430 | 0.1794 | 0.2910 | 0.0001 |
| sweep 0.7, size-matched random | TCR-IMAP vs pair-frequency prior | canonical CDRs only | 174 | 0.5996 | 0.3566 | 0.2430 | 0.1794 | 0.2911 | 0.0001 |
| sweep 0.7, similarity-filtered | TCR-IMAP vs pair-frequency prior | all rows (reported) | 174 | 0.3271 | 0.2997 | 0.0274 | 0.0011 | 0.0542 | 0.0414 |
| sweep 0.7, similarity-filtered | TCR-IMAP vs pair-frequency prior | canonical CDRs only | 174 | 0.3273 | 0.2998 | 0.0275 | 0.0011 | 0.0542 | 0.0412 |

**What this covers.** 14 comparisons could be rescored: TCR-IMAP against the pair-frequency prior on the
pairwise task under all five protocols, TCR-IMAP against TEIM-Res under tcr70, TCR-IMAP against Struct-TCS and against
Seq-TCS on the per-residue task over the 88-complex comparator subset, and the three sweep thresholds in both the
similarity-filtered and the size-matched-random regime. Across those 14 the largest change in the paired
difference is **−0.0014** (pairwise, tcr70, TCR-IMAP vs TEIM-Res); every other comparison moves by at most
0.0002, and no interval crosses zero that did not already. TEIM-Res’s own macro-AUPRC under tcr70 moves from
0.2725 to 0.2740. Within this listed set the FR4 rows change no verdict,
which is a statement about these 14 comparisons and not about every result in the paper.

**What this does not cover.** The **predeclared primary comparison**, per-residue macro-AUPRC under tcr70, TCR-IMAP
versus Struct-TCS, over the full 174-complex benchmark, **cannot be rescored at all**. That runner stores per-complex
metrics and no prediction array, and writes no checkpoint, so there is nothing to restrict to canonical CDR rows; the
files under `results/checkpoints/` are from the superseded 28-position candidate set under partitions that are not the
frozen protocols, and are not the reported runs (`results/review/FR4_CHECKPOINT_AUDIT.md`). Completing this one
comparison would require retraining the per-residue benchmark with predictions persisted. It is reported as
outstanding, **not** estimated from the rows above: the 88-complex comparator-subset figures are a different
evaluation on a different set of complexes and are not a proxy for it. The gap is in one robustness check; the primary
comparison itself was computed over its full predeclared candidate universe.
