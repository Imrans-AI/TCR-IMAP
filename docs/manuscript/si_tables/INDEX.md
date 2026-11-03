# Supporting Information tables: index

Built by `build_si_tables.py` (this directory) from committed tables, reports and stored per-complex arrays. No
training, inference, scoring or GPU. Recomputed statistics use the existing pooled clu70 cluster bootstrap
(`src/viz/figstats.pooled_cluster_bootstrap`, 20,000 resamples, seed 42) on stored per-complex arrays and are asserted
equal to the committed values; the build stops on any mismatch. Log: `build.log`. Numbering follows the task order
(S1–S8), not yet citation order; these tables are distinct from Figures S1–S9.

| Table | Title | Caption (short) | Source files | CSV |
|---|---|---|---|---|
| [S1](S1_benchmark.md) | Full benchmark: every method under every protocol | Per-complex macro-AUPRC and ROC AUC of every method under every protocol, n complexes and clusters, paired differences from TCR-IMAP (per-residue, pairwise, comparator subset, AngleNet) | `percomplex_{protocol}_p84v2.npz`, `benchmark_p84v2.csv`, `protocol_main_p84v2.csv`, `pairwise_ours_p84v2.npz`, `teim_res_pairwise_p84v2.npz`, `structcs_setting_p84v2_*.csv`, `angle_runs_p84v2.jsonl`; `figures_p84v2/stats_20k_seed42.tsv`, `validation_fix/F_structcs_differences_20k.tsv` | `S1_benchmark_anglenet.csv`, `S1_benchmark_comparator_subset.csv`, `S1_benchmark_per_residue_pairwise.csv` |
| [S2](S2_ablations.md) | Ablations of TCR-IMAP (tcr70, pairwise task) | All eight ablations: Δ, 95 % CI, Bonferroni-simultaneous 99.375 % interval, raw and Holm-adjusted p, survival | `_partial/p84v2/29b_ablations/*.npz`, `ablations_p84v2.csv`; `figures_p84v2/stats_20k_seed42.tsv` | `S2_ablations.csv` |
| [S3](S3_hyperparameters.md) | Architectures, hyperparameters and parameter counts | Layers, widths, heads, dropout, optimizer, learning rate, weight decay, epochs, early stopping, seeds and parameter counts of TCR-IMAP (both variants), Struct-TCS, Seq-TCS, TEIM-Res and AngleNet | `scripts/19`, `21`, `24b`, `28b`, `53_benchmark_p84v2`, `20`, `36b`, `16b`; parameter counts by instantiation | `S3_hyperparameters.csv` |
| [S4](S4_protocol_folds.md) | Partitioning protocols: fold composition | Per protocol and fold: n training, n test, clusters on each side, maximum and mean train/test similarity, duplicate families split; allele/intersection coverage limitation | `protocols.json`, `complex_metadata.tsv`, `tcr_identity_matrix.npy`, `leakage_audit_dual.tsv`, `percomplex_maxsim.tsv`, `validation_fix/B2_split_dryrun.log` | `S4_protocol_folds.csv` |
| [S5](S5_identity_sweep.md) | Identity-threshold sweep: composition and results | Per threshold: n test, n training per fold, achieved maximum and mean similarity, the `mean_max_sim` upper bound, Δ, CI, p | `similarity_curve_p84v2.csv`, `similarity_curve_composition{,_folds}_p84v2.csv` | `S5_identity_sweep.csv` |
| [S6](S6_md_selection_qc.md) | Molecular dynamics: system selection and per-replicate quality control | The four systems (selection basis, allele, peptide, resolution, p84v2 rank) and, per replicate, steps, frames, LINCS warnings, temperature, pressure, density, backbone RMSD, TCR–pMHC contact retention | `md/MD_SYSTEM_SELECTION_CRITERIA.md`, `percomplex_tcr70_paired_p84v2.tsv`, `characterisation_fix/S2_qc_all12.log`, `md_characterisation_*.csv`, `md_occupancy_docked_frames.csv` | `S6_md_replicates.csv`, `S6_md_systems.csv` |
| [S7](S7_md_occupancy.md) | Molecular dynamics: occupancy statistics | Pooled, per-system, per-stratum and per-replicate Spearman ρ for TCR-IMAP and the prior with counts; top-k enrichment; crystal-positive occupancy; docked-frame analysis; score–score correlation | `md_occupancy_{rho,strata,topk,check1_positive,docked_rho}.csv`, `md_score_score_spearman.csv` | `S7_md_check1.csv`, `S7_md_docked_rho.csv`, `S7_md_rho.csv`, `S7_md_score_score.csv`, `S7_md_strata.csv`, `S7_md_topk.csv` |
| [S8](S8_dataset_descriptors.md) | Dataset descriptors | Exact twins, twin pairs differing in peptide, deduplication (132 complexes, median NN 0.967), helix share 88.3 %, β-sheet floor, each with its computation | `tcr_identity_matrix.npy`, `complex_metadata.tsv`, `identity_bound/handoff6_statements.tsv`, `gdomain182_region_loop_summary_p84v2.csv` | `S8_dataset_descriptors.csv` |

## Pending

- none
