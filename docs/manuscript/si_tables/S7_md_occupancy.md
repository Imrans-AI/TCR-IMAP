# Table S7. Molecular dynamics: occupancy statistics

**Table S7.** Spearman correlation between score and solution-phase contact occupancy among crystal-negative residue
pairs, for TCR-IMAP and the pair-frequency prior (the predeclared control), with the occupancy counts each correlation
rests on. (A) Pooled and per system, overall and by crystal-frame distance stratum. (B) Per replicate. (C) Top-k
enrichment. (D) Occupancy of crystal-positive pairs (check 1). (E) Docked-frame secondary analysis. (F) Correlation
between the two scores.

## A. Primary statistic, pooled and per system, by stratum

| System | Stratum | n pairs | n non-zero | n ≥ 0.25 | ρ TCR-IMAP | p (approx.) | ρ prior | p (approx.) |
|---|---|---|---|---|---|---|---|---|
| **pooled (all four)** | all | 19,957 | 381 | 34 | 0.216 | < 0.001 | 0.441 | < 0.001 |
|  | 4.5–6 Å | 120 | 110 | 25 | 0.241 | 0.008 | 0.179 | 0.046 |
|  | 6–8 Å | 423 | 177 | 9 | 0.303 | < 0.001 | 0.309 | < 0.001 |
|  | > 8 Å | 19,414 | 94 | 0 | 0.103 | < 0.001 | 0.208 | < 0.001 |
| 5yxn | all | 4,920 | 119 | 6 | 0.236 | < 0.001 | 0.471 | < 0.001 |
|  | 4.5–6 Å | 38 | 33 | 4 | 0.301 | 0.066 | 0.284 | 0.083 |
|  | 6–8 Å | 118 | 50 | 2 | 0.304 | 0.001 | 0.284 | 0.002 |
|  | > 8 Å | 4,764 | 36 | 0 | 0.121 | < 0.001 | 0.254 | < 0.001 |
| 7n1f | all | 4,598 | 83 | 8 | 0.219 | < 0.001 | 0.488 | < 0.001 |
|  | 4.5–6 Å | 24 | 21 | 7 | 0.305 | 0.143 | 0.071 | 0.739 |
|  | 6–8 Å | 92 | 41 | 1 | 0.323 | 0.002 | 0.305 | 0.003 |
|  | > 8 Å | 4,482 | 21 | 0 | 0.109 | < 0.001 | 0.265 | < 0.001 |
| 3vxs | all | 5,010 | 73 | 6 | 0.195 | < 0.001 | 0.439 | < 0.001 |
|  | 4.5–6 Å | 22 | 22 | 5 | 0.336 | 0.122 | 0.312 | 0.153 |
|  | 6–8 Å | 104 | 39 | 1 | 0.448 | < 0.001 | 0.460 | < 0.001 |
|  | > 8 Å | 4,884 | 12 | 0 | 0.072 | < 0.001 | 0.151 | < 0.001 |
| 8eo8 | all | 5,429 | 106 | 14 | 0.213 | < 0.001 | 0.393 | < 0.001 |
|  | 4.5–6 Å | 36 | 34 | 9 | 0.151 | 0.386 | 0.023 | 0.892 |
|  | 6–8 Å | 109 | 47 | 5 | 0.176 | 0.065 | 0.219 | 0.023 |
|  | > 8 Å | 5,284 | 25 | 0 | 0.110 | < 0.001 | 0.183 | < 0.001 |

## B. Per replicate (all crystal-negative pairs)

| System | Rep | n pairs | n non-zero | n ≥ 0.25 | ρ TCR-IMAP | ρ prior |
|---|---|---|---|---|---|---|
| 5yxn | 1 | 4,920 | 84 | 2 | 0.198 | 0.386 |
| 5yxn | 2 | 4,920 | 66 | 10 | 0.185 | 0.401 |
| 5yxn | 3 | 4,920 | 49 | 5 | 0.155 | 0.321 |
| 7n1f | 1 | 4,598 | 49 | 5 | 0.167 | 0.408 |
| 7n1f | 2 | 4,598 | 64 | 8 | 0.194 | 0.435 |
| 7n1f | 3 | 4,598 | 47 | 6 | 0.169 | 0.377 |
| 3vxs | 1 | 5,010 | 50 | 6 | 0.164 | 0.380 |
| 3vxs | 2 | 5,010 | 48 | 2 | 0.162 | 0.360 |
| 3vxs | 3 | 5,010 | 51 | 9 | 0.165 | 0.383 |
| 8eo8 | 1 | 5,429 | 65 | 14 | 0.166 | 0.324 |
| 8eo8 | 2 | 5,429 | 85 | 12 | 0.190 | 0.360 |
| 8eo8 | 3 | 5,429 | 79 | 15 | 0.180 | 0.332 |

## C. Top-k enrichment (k = the system's number of crystal-positive pairs)

| System | Score | k | Top-k pairs with occupancy ≥ 0.25 (fraction) | Base rate (all crystal-negative pairs) | Enrichment | Tie in score at rank k |
|---|---|---|---|---|---|---|
| 5yxn | TCR-IMAP | 36 | 3 / 36 (0.083) | 0.0012 | 68× | no |
|  | pair-frequency prior | 36 | 3 / 36 (0.083) | 0.0012 | 68× | no |
| 7n1f | TCR-IMAP | 22 | 3 / 22 (0.136) | 0.0017 | 78× | no |
|  | pair-frequency prior | 22 | 2 / 22 (0.091) | 0.0017 | 52× | yes |
| 3vxs | TCR-IMAP | 30 | 1 / 30 (0.033) | 0.0012 | 28× | no |
|  | pair-frequency prior | 30 | 1 / 30 (0.033) | 0.0012 | 28× | no |
| 8eo8 | TCR-IMAP | 31 | 5 / 31 (0.161) | 0.0026 | 63× | no |
|  | pair-frequency prior | 31 | 7 / 31 (0.226) | 0.0026 | 88× | yes |

## D. Crystal-positive occupancy (check 1)

| System | Rep | n crystal-positive pairs | Mean | Median | Min | Occupancy ≥ 0.25 |
|---|---|---|---|---|---|---|
| 5yxn | 1 | 36 | 0.168 | 0.065 | 0.000 | 6 / 36 |
| 5yxn | 2 | 36 | 0.431 | 0.490 | 0.003 | 24 / 36 |
| 5yxn | 3 | 36 | 0.393 | 0.334 | 0.000 | 20 / 36 |
| 5yxn | system | 36 | 0.331 | 0.312 | 0.002 | 20 / 36 |
| 7n1f | 1 | 22 | 0.321 | 0.203 | 0.000 | 8 / 22 |
| 7n1f | 2 | 22 | 0.525 | 0.495 | 0.007 | 16 / 22 |
| 7n1f | 3 | 22 | 0.594 | 0.699 | 0.009 | 17 / 22 |
| 7n1f | system | 22 | 0.480 | 0.531 | 0.014 | 16 / 22 |
| 3vxs | 1 | 30 | 0.586 | 0.788 | 0.000 | 20 / 30 |
| 3vxs | 2 | 30 | 0.634 | 0.847 | 0.003 | 22 / 30 |
| 3vxs | 3 | 30 | 0.573 | 0.625 | 0.002 | 22 / 30 |
| 3vxs | system | 30 | 0.598 | 0.674 | 0.006 | 24 / 30 |
| 8eo8 | 1 | 31 | 0.519 | 0.537 | 0.003 | 20 / 31 |
| 8eo8 | 2 | 31 | 0.463 | 0.552 | 0.001 | 17 / 31 |
| 8eo8 | 3 | 31 | 0.458 | 0.441 | 0.000 | 19 / 31 |
| 8eo8 | system | 31 | 0.480 | 0.598 | 0.008 | 20 / 31 |

## E. Docked-frame secondary analysis (system occupancy over docked frames of contributing replicates)

| Reading | System | Replicates contributing | n pairs | n non-zero | ρ TCR-IMAP | ρ prior | Interpretable |
|---|---|---|---|---|---|---|---|
| count | 5yxn | 1 | 4,920 | 55 | 0.170 | 0.390 | yes |
| count | 7n1f | 3 | 4,598 | 83 | 0.219 | 0.488 | yes |
| count | 3vxs | 3 | 5,010 | 73 | 0.195 | 0.439 | yes |
| count | 8eo8 | 3 | 5,429 | 106 | 0.213 | 0.393 | yes |
| count | 4 contributing systems | n/a | 19,957 | 317 | 0.200 | 0.424 | yes |
| identity | 5yxn | 0 | 4,920 | n/a | n/a | n/a | no (no contributing replicate) |
| identity | 7n1f | 2 | 4,598 | 50 | 0.172 | 0.402 | yes |
| identity | 3vxs | 3 | 5,010 | 59 | 0.176 | 0.407 | yes |
| identity | 8eo8 | 0 | 5,429 | n/a | n/a | n/a | no (no contributing replicate) |
| identity | 2 contributing systems | n/a | 9,608 | 109 | 0.174 | 0.405 | yes |

## F. Score–score Spearman correlation among crystal-negative pairs

| System | n pairs | ρ (TCR-IMAP vs prior) |
|---|---|---|
| 5yxn | 4,920 | 0.341 |
| 7n1f | 4,598 | 0.334 |
| 3vxs | 5,010 | 0.349 |
| 8eo8 | 5,429 | 0.395 |
| all four | 19,957 | 0.355 |

**Notes.**
- Occupancy of a pair: fraction of the 9,001 frames between 20 and 200 ns in which its minimum heavy-atom distance is
  ≤ 4.5 Å (the crystal-label definition); system occupancy is the mean of the three replicates. Pair universe: resolved
  pairs of the p84v2 label array (75 TCR rows × 84 MHC columns); crystal-negative = not a contact in the crystal
  structure. Scores: TCR-IMAP's out-of-fold tcr70 pairwise predictions and the within-fold pair-frequency prior.
- p: within-system permutation test, 10,000 permutations of the scores among a system's crystal-negative pairs (within
  system for the pooled rows), seed 42, two-sided, **approximate** because pairs within a system are not independent;
  p is floored at 1/10,001 and shown as < 0.001. Per-replicate rows carry no p. No test of the difference between the
  two scores was made.
- Strata: the pair's minimum heavy-atom distance in the crystal structure, 4.5–6, 6–8 and > 8 Å (fixed before the run).
- Every row in A and B rests on at least 10 non-zero occupancies and is flagged `interpretable`; in E, the identity
  reading has no contributing replicate for 5yxn and 8eo8 (no rho; reported, not interpreted).
- C: enrichment = fraction of top-k pairs with occupancy ≥ 0.25 ÷ base rate; it rests on 1–7 pairs per system. A tie in
  score at rank k is flagged; the top-k set uses a stable descending sort in row-major pair order.
- D: no numeric acceptance criterion was predeclared for check 1.
- E: docked frame, count reading (predeclared, primary for this analysis): at least half as many TCR–pMHC heavy-atom
  contacts as at t = 0; identity reading: at least half of the specific t = 0 contact pairs present. A replicate
  contributes only with ≥ 1,000 docked frames.

**Sources.** `results/tables/md_occupancy_rho.csv` (A overall rows, B), `md_occupancy_strata.csv` (A strata),
`md_occupancy_topk.csv` (C), `md_occupancy_check1_positive.csv` (D), `md_occupancy_docked_rho.csv` (E),
`md_score_score_spearman.csv` (F); all written once by the single authorised run (manifest `md_occupancy_manifest.csv`)
and read with `float_precision="round_trip"`. Definitions: `md/MD_SYSTEM_SELECTION_CRITERIA.md` (analysis definitions of
20 Oct; revisions and decisions of 24 Oct).
