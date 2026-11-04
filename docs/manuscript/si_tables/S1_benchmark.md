# Table S1. Full benchmark: every method under every protocol

**Table S1.** Per-complex macro-averaged AUPRC and ROC AUC of every method under every partitioning protocol, with the
number of test complexes and of TCR clusters (clu70) they span, and the paired difference from TCR-IMAP. (A) Per-residue
task. (B) Pairwise task. (C) The 88 structures shared with the comparator's dataset. (D) Docking-angle regression
(AngleNet, a separate regressor; not TCR-IMAP).

## A. Per-residue task (cluster-disjoint validation; the reported configuration)

| Protocol | n complexes | n clusters | Method | Macro-AUPRC | Macro-AUC | TCR-IMAP − method [95 % CI] | p |
|---|---|---|---|---|---|---|---|
| random | 174 | 69 | TCR-IMAP | 0.8061 | 0.9309 | n/a | n/a |
|  |  |  | Struct-TCS | 0.8092 | 0.9280 | −0.0030 [−0.0227, +0.0207] | 0.798 |
|  |  |  | Seq-TCS | 0.7736 | 0.9074 | +0.0325 [+0.0132, +0.0559] | 0.001 |
|  |  |  | Position-frequency prior | 0.7457 | 0.9086 | +0.0604 [+0.0364, +0.0837] | < 0.001 |
|  |  |  | Nearest-neighbor transfer | 0.7281 | 0.8692 | +0.0780 [+0.0398, +0.1252] | < 0.001 |
| intersection (legacy) | 30 | 22 | TCR-IMAP | 0.7996 | 0.9293 | n/a | n/a |
|  |  |  | Struct-TCS | 0.7300 | 0.8944 | +0.0696 [+0.0148, +0.1161] | 0.015 |
|  |  |  | Seq-TCS | 0.7079 | 0.8797 | +0.0917 [+0.0394, +0.1436] | < 0.001 |
|  |  |  | Position-frequency prior | 0.7213 | 0.8988 | +0.0783 [+0.0431, +0.1127] | < 0.001 |
|  |  |  | Nearest-neighbor transfer | 0.5704 | 0.7796 | +0.2292 [+0.1527, +0.2943] | < 0.001 |
| allele | 64 | 38 | TCR-IMAP | 0.7764 | 0.9237 | n/a | n/a |
|  |  |  | Struct-TCS | 0.7258 | 0.8875 | +0.0507 [+0.0187, +0.0822] | 0.001 |
|  |  |  | Seq-TCS | 0.6734 | 0.8615 | +0.1030 [+0.0716, +0.1382] | < 0.001 |
|  |  |  | Position-frequency prior | 0.7467 | 0.9096 | +0.0298 [−0.0163, +0.0775] | 0.205 |
|  |  |  | Nearest-neighbor transfer | 0.5848 | 0.7840 | +0.1916 [+0.1001, +0.2820] | < 0.001 |
| temporal | 60 | 32 | TCR-IMAP | 0.6936 | 0.8861 | n/a | n/a |
|  |  |  | Struct-TCS | 0.6907 | 0.8789 | +0.0030 [−0.0410, +0.0602] | 0.874 |
|  |  |  | Seq-TCS | 0.6474 | 0.8485 | +0.0462 [+0.0043, +0.1064] | 0.028 |
|  |  |  | Position-frequency prior | 0.6830 | 0.8841 | +0.0107 [−0.0257, +0.0515] | 0.584 |
|  |  |  | Nearest-neighbor transfer | 0.4514 | 0.7112 | +0.2423 [+0.1925, +0.2913] | < 0.001 |
| tcr70 | 174 | 69 | TCR-IMAP | 0.7281 | 0.8968 | n/a | n/a |
|  |  |  | Struct-TCS | 0.6969 | 0.8740 | +0.0312 [−0.0049, +0.0650] | 0.088 |
|  |  |  | Seq-TCS | 0.6481 | 0.8453 | +0.0799 [+0.0455, +0.1161] | < 0.001 |
|  |  |  | Position-frequency prior | 0.7196 | 0.8982 | +0.0085 [−0.0162, +0.0371] | 0.546 |
|  |  |  | Nearest-neighbor transfer | 0.4596 | 0.7182 | +0.2685 [+0.2368, +0.3004] | < 0.001 |

## B. Pairwise task

| Protocol | n complexes | n clusters | Method | Macro-AUPRC | Macro-AUC | TCR-IMAP − method [95 % CI] | p | TCR-IMAP higher |
|---|---|---|---|---|---|---|---|---|
| random | 174 | 69 | TCR-IMAP | 0.6230 | 0.9912 | n/a | n/a | n/a |
|  |  |  | Pair-frequency prior | 0.3559 | 0.9741 | +0.2671 [+0.2069, +0.3126] | < 0.001 | 151/174 |
| intersection (legacy) | 30 | 22 | TCR-IMAP | 0.4406 | 0.9876 | n/a | n/a | n/a |
|  |  |  | Pair-frequency prior | 0.3182 | 0.9730 | +0.1223 [+0.0508, +0.2021] | < 0.001 | 22/30 |
| allele | 64 | 38 | TCR-IMAP | 0.4721 | 0.9832 | n/a | n/a | n/a |
|  |  |  | Pair-frequency prior | 0.3306 | 0.9689 | +0.1415 [+0.0577, +0.2241] | < 0.001 | 48/64 |
| temporal | 60 | 32 | TCR-IMAP | 0.3200 | 0.9768 | n/a | n/a | n/a |
|  |  |  | Pair-frequency prior | 0.2892 | 0.9354 | +0.0308 [+0.0012, +0.0663] | 0.041 | 40/60 |
| tcr70 | 174 | 69 | TCR-IMAP | 0.3598 | 0.9809 | n/a | n/a | n/a |
|  |  |  | Pair-frequency prior | 0.3090 | 0.9601 | +0.0507 [+0.0215, +0.0841] | < 0.001 | 112/174 |
|  |  |  | TEIM-Res | 0.2725 | 0.8205 | +0.0873 [+0.0592, +0.1149] | < 0.001 | 143/174 |

## C. Comparator subset: 88 structures shared with the Struct-TCS dataset (per-residue; cluster-disjoint validation)

| Cutoff | Setting | n complexes | n clusters | Comparator | TCR-IMAP macro-AUPRC | Comparator macro-AUPRC | TCR-IMAP − comparator [95 % CI] | p |
|---|---|---|---|---|---|---|---|---|
| 4.5 Å | single random split | 22 | 17 | Struct-TCS | 0.8108 | 0.7839 | +0.0269 [−0.0102, +0.0707] | 0.171 |
| 4.5 Å | single random split | 22 | 17 | Seq-TCS | 0.8108 | 0.7375 | +0.0733 [+0.0308, +0.1155] | < 0.001 |
| 4.5 Å | tcr70 | 88 | 40 | Struct-TCS | 0.7297 | 0.6241 | +0.1057 [+0.0705, +0.1357] | < 0.001 |
| 4.5 Å | tcr70 | 88 | 40 | Seq-TCS | 0.7297 | 0.6023 | +0.1274 [+0.0884, +0.1618] | < 0.001 |
| 6 Å | single random split | 22 | 17 | Struct-TCS | 0.7621 | 0.8570 | −0.0949 [−0.1335, −0.0517] | < 0.001 |
| 6 Å | single random split | 22 | 17 | Seq-TCS | 0.7621 | 0.8333 | −0.0712 [−0.1195, −0.0215] | 0.004 |
| 6 Å | tcr70 | 88 | 40 | Struct-TCS | 0.8166 | 0.7202 | +0.0964 [+0.0684, +0.1235] | < 0.001 |
| 6 Å | tcr70 | 88 | 40 | Seq-TCS | 0.8166 | 0.7392 | +0.0774 [+0.0545, +0.1045] | < 0.001 |

## D. Docking-angle regression (AngleNet; one row per partition and training seed)

| Partition | Seed | n | MAE (°) | Median-predictor MAE (°) | Pearson r | Spearman ρ | R² |
|---|---|---|---|---|---|---|---|
| random | 42 | 172 | 5.40 | 11.00 | 0.776 | 0.775 | +0.556 |
| random | 43 | 172 | 5.80 | 11.00 | 0.732 | 0.746 | +0.454 |
| random | 44 | 172 | 6.08 | 11.00 | 0.713 | 0.709 | +0.420 |
| allele | 42 | 172 | 11.18 | 11.00 | 0.133 | 0.068 | −0.382 |
| allele | 43 | 172 | 11.22 | 11.00 | 0.196 | 0.093 | −0.310 |
| allele | 44 | 172 | 12.01 | 11.00 | 0.159 | 0.108 | −0.216 |
| trbv | 42 | 172 | 12.50 | 11.00 | 0.202 | 0.205 | −0.394 |
| trbv | 43 | 172 | 13.08 | 11.00 | 0.110 | 0.083 | −0.468 |
| trbv | 44 | 172 | 13.43 | 11.00 | 0.128 | 0.140 | −0.407 |

**Notes.**
- Macro-AUPRC and macro-AUC are means over test complexes of the per-complex value; a complex with no positive or no
  negative resolved candidate has no defined AUPRC or AUC and is omitted from that mean (none is omitted in any cell of A).
- Differences are paired: TCR-IMAP − method, per complex, then averaged. 95 % CI and two-sided p from the pooled clu70
  cluster bootstrap (`src/viz/figstats.pooled_cluster_bootstrap`), **20,000 resamples, seed 42, recomputed by
  `build_si_tables.py` from the stored per-complex arrays** and asserted equal to the committed values; p is floored at
  1/20,000 and reported as < 0.001 at or below 0.0005 (project convention). The p values in the
  `vs_*` columns of `benchmark_p84v2.csv` and `structcs_setting_p84v2_*.csv` are the run's own 2,000-resample values and
  differ slightly; the 20,000-resample values here are the ones reported in the paper.
- "TCR-IMAP higher": complexes on which TCR-IMAP's AUPRC exceeds the comparator's.
- n clusters: distinct clu70 clusters among the evaluated complexes (the bootstrap unit).
- Pairwise ROC AUC is computed here by scoring the stored prediction arrays in
  `results/tables/pairmaps_{protocol}_p84v2.npz`, exactly as the macro-AUPRC in the same rows is. Those arrays were
  verified against the frozen labels, masks, identifiers and cluster assignments before use, and the recomputed
  tcr70 values reproduce the committed `pairwise_ours_p84v2.npz` AUC and the macro-AUPRC of every protocol.
- D: R² is relative to the mean predictor; the MAE baseline is the median predictor. AngleNet uses its own five-fold
  random, allele-holdout and TRBV-holdout partitions (not the frozen protocols); 172 of 174 complexes (3vxu and 3w0w have
  no docking angle).

**Sources.** A: `results/tables/percomplex_{protocol}_p84v2.npz` (SHA-256 verified against
`results/review/validation_fix/clusterdisjoint_run1/SHA256SUMS`), means checked against `results/tables/benchmark_p84v2.csv`,
differences checked against `results/review/figures_p84v2/stats_20k_seed42.tsv`. B: `results/tables/protocol_main_p84v2.csv`
(means), `stats_20k_seed42.tsv` (differences, wins; recomputed for tcr70 from `percomplex_tcr70_paired_p84v2.tsv`),
`pairwise_ours_p84v2.npz` and `teim_res_pairwise_p84v2.npz` (tcr70 AUC, TEIM-Res). The pairwise differences under
random, intersection (legacy), allele and temporal are taken from `stats_20k_seed42.tsv` as committed (their per-complex arrays were
derived in memory by `scripts/55_stats_p84v2.py` and not persisted). C: `results/tables/structcs_setting_p84v2_{4p5A,6A}.csv`,
`results/tables/_partial/p84v2/31b_structcs_{45A,6A}/*.npz`, checked against
`results/review/validation_fix/F_structcs_differences_20k.tsv` (cluster-disjoint rows). D: `results/tables/angle_runs_p84v2.jsonl`.
