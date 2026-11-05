# Figure manifest

**Final figures**, copied into `figures/main/` and `figures/supplementary/` on 26 October 2026 (branch
`figure-audit`). The originals in `figures/` are untouched, and every script still writes to its original path. **The
copies are not regenerated automatically:** if a figure is rebuilt, re-copy it here.

Each copy is byte-identical to its source at the time of copying. The source is the state of `figS-md-v3` @ `b82cb3d`,
which is `master` plus the author's six SI-figure commits.

**Re-copied 26 October 2026 (branch `main-figure-fixes`):** Figures 2–7 were rebuilt for six presentation fixes
(FIGURE_AUDIT_REPORT.md, "Addendum 2") and re-copied; each copy is again byte-identical to its source. Figure 2 was
rebuilt only for spelling (neighbour → neighbor); its inputs and every value are unchanged. The SI copies were not
rebuilt.

## Main figures (`figures/main/`)

| New name | Source file (`figures/`) | Script | Tables read | Description |
|---|---|---|---|---|
| `Figure_1.png/.pdf/.svg` | n/a (author-approved asset, not a script output; see `figures/main/Figure_1.SOURCE.md`) | `figures/source/fig1_arial_20261102/design_figure.py` (editable source, Arial revision of 2 Nov 2026; the approved DejaVu Sans source is preserved at `figures/source/fig1_approved_20261102/`; `render_structure.py` for the PyMOL panel) | `results/tables/pairmaps_tcr70_p84v2.npz`, `data/processed/complex_metadata.tsv`, `data/structures/stcrdab/imgt/1ao7.pdb` | Structural benchmark and sequence-based contact prediction. (A) benchmark, contact labels, aligned sequences and the controlled evaluation, with the 1AO7 complex; (B) the pairwise TCR-IMAP architecture, with 1AO7's crystal labels and its stored out-of-fold tcr70 prediction. **Installed 2 Nov 2026, font-only revision to Arial the same day; do not regenerate.** `scripts/70_fig1_workflow.py` is superseded and now refuses to run without `--allow-superseded`; the manuscript builder verifies the figure's SHA-256 and aborts on a mismatch |
| `TOC_graphic.png/.tif/.pdf/.svg` | Author-approved exports, installed 5 November 2026; see `figures/main/TOC_graphic.SOURCE.md` | `figures/source/toc_reference_approved_20261105/build_toc.py` (editable composition; installed assets are not regenerated automatically) | `results/tables/pairmaps_tcr70_p84v2.npz`, `results/tables/similarity_curve_p84v2_authoritative.csv`; molecular rendering from the corresponding experimental complex | TOC graphic, 8.25 × 4.45 cm; PNG/TIFF at 900 dpi. Three stages: molecular context and sequence fragments, stored contact prediction, measured similarity-threshold evaluation. Arial 6 pt, Courier New sequence tokens; only stage labels bold; no legends. Supersedes the 4 November TOC and its old scripts. Author-approved files copied byte-for-byte; previous version preserved. |
| `Figure_2.png/.pdf` | `fig2_dataset_p82v1_v1` | `scripts/38b_fig2_dataset_p82v1.py` | `data/processed/complex_metadata.tsv`, `results/tables/percomplex_maxsim.tsv`, `results/tables/leakage_audit_dual.tsv`, `results/tables/percomplex_{protocol}_p82v1.npz` | Dataset redundancy and leakage. **Keeps its p82v1 name by design:** its inputs are the frozen metadata and splits, and it is byte-identical under p84v2 (the project notes; FIGURES_P84V2_REPORT.md) |
| `Figure_3.png/.pdf` | `fig3_inflation_p84v2_v1` | `scripts/55_fig3_inflation_p84v2.py` | `results/tables/similarity_curve_p84v2.csv`, `angle_runs_p84v2.jsonl` | Inflation with sequence similarity (identity sweep) and the AngleNet docking-angle collapse under grouped evaluation |
| `Figure_4.png/.pdf` | `fig4_benchmark_p84v2_v1` | `scripts/55_fig4_benchmark_p84v2.py` | `benchmark_p84v2.csv`, `protocol_main_p84v2.csv`, `percomplex_{protocol}_p84v2.npz`, `percomplex_tcr70_paired_p84v2.tsv`, `pairwise_ours_p84v2.npz`, `teim_res_pairwise_p84v2.npz`, `cluster_macro_p84v2.csv`, `structcs_setting_p84v2_{4p5A,6A}.csv`, `results/review/figures_p84v2/stats_20k_seed42.tsv`, `complex_metadata.tsv` | Benchmark against published comparators: per-residue by protocol, paired differences, pairwise vs prior and TEIM-Res, cluster-macro sensitivity, comparator subset |
| `Figure_5.png/.pdf` | `fig5_footprint_p84v2_v1` | `scripts/55_fig5_footprint_p84v2.py`, with PyMOL panels from `scripts/55_pymol_fig5_p84v2.py` (`figures/pymol/p84v2_panel{A_complex,B_footprint}.png`) | `contact_freq_gdomain182_descriptive_p84v2.npz`, `gdomain182_region_classification_p84v2.csv`, `germline_diagonal_p84v2.csv`; PyMOL: `stcrdab_index.tsv`, the residue manifest, IMGT structures | Observed structural contact footprint on the α1/α2 window: (A) contact frequency, 75 TCR rows × 84 window columns; (B) the 1ao7 MHC platform by relative contact frequency (PyMOL); (C) the 1ao7 TCR–pMHC complex with CDR loops (PyMOL; intentional, loop selections verified against the residue manifest, `fig5_complex_cdr_selection_check.tsv`); (D) contact rate by loop × MHC region; (E) per-cluster α1 − α2 difference |
| `Figure_6.png/.pdf` | `fig6_ablation_p84v2_v1` | `scripts/55_fig6_ablation_p84v2.py`, with PyMOL panels from `scripts/55_pymol_fig6_p84v2.py` (`figures/pymol/p84v2_fig6_{3qfj,3qeq}_{observed,predicted}.png`) | `ablations_p84v2.csv`, ablation `{key}.npz` and `full.npz`, `pairwise_ours_p84v2.npz`, `percomplex_tcr70_paired_p84v2.tsv`, `stats_20k_seed42.tsv`, `figures/pymol/fig6_sites_p84v2.tsv`, `complex_metadata.tsv` | Ablations (Holm-adjusted) and the two showcase complexes, 3qfj and 3qeq |
| `Figure_7.png/.pdf` | `fig7_md_occupancy_v1` | `scripts/56_fig7_md_occupancy.py` | the 8 `md_occupancy_*` tables, `md_characterisation_{summary,timeseries,manifest}.csv`, `pairmaps_tcr70_p84v2.npz` (sha256-checked against the run manifest), `data/processed/pair_maps_p84v2.npz` | MD occupancy of crystal-negative pairs vs TCR-IMAP and the prior: (A) ρ by stratum; (B) score vs occupancy; (C) top-k; (D) check 1 with docked fractions; (E) backbone RMSD |

## Supplementary figures (`figures/supplementary/`)

**Numbering:** S1–S3 follow the order the audit task lists them, and S4–S9 the order of the SI-figures task, not
manuscript citation order. Renumber when the SI order is fixed (after the Results text exists).

| New name | Source file (`figures/`) | Script | Tables read | Description |
|---|---|---|---|---|
| `Figure_S1.png/.pdf` | `figS_peptide_p84v2_v1` | `scripts/55_figS_peptide_p84v2.py` | `results/tables/peptide_aware_p84v2.csv`, `stats_20k_seed42.tsv` | Effect of adding peptide sequence to TCR-IMAP: a null in every scope |
| `Figure_S2.png/.pdf` | `figS_benchmark_distributions_p84v2_v1` | `scripts/55_figS_benchmark_dist_p84v2.py` | `benchmark_p84v2.csv`, `percomplex_{protocol}_p84v2.npz` | Per-complex per-residue AUPRC distributions under each partitioning protocol |
| `Figure_S3.png/.pdf` | `figS_md_characterisation_v1` | `scripts/56_figS_md_characterisation.py` | `md_characterisation_{summary,timeseries,manifest}.csv`, `md_rmsf_per_residue.csv` | MD interface characterisation: interface RMSD, buried SASA, per-loop RMSF, H-bonds (descriptive) |
| `Figure_S4.png/.pdf` | `figS_validation_split_p84v2_v1` | `scripts/57_figS_validation_split.py` | `validation_fix/randomsplit_deterministic/percomplex_{protocol}_p84v2.npz`, `percomplex_{protocol}_p84v2.npz` (both sha256-checked), `validation_fix/B3_macro_auprc.tsv`, `B3_key_differences_20k.tsv`, `queue_console_randomsplit_deterministic.log`, `B2_split_dryrun.log` | Validation-split sensitivity: per-residue macro-AUPRC under the random vs cluster-disjoint split, tcr70 key differences under both, clusters straddling model selection per fold |
| `Figure_S5.png/.pdf` | `figS_cluster_macro_p84v2_v1` | `scripts/57_figS_cluster_macro.py` | `percomplex_tcr70_p84v2.npz`, `percomplex_tcr70_paired_p84v2.tsv`, `pairwise_ours_p84v2.npz`, `figures_p84v2/stats_20k_seed42.tsv`, `cluster_macro_p84v2.csv` | Cluster-macro sensitivity of the three primary tcr70 comparisons (complex-macro primary) |
| `Figure_S6.png/.pdf` | `figS_comparator_subset_seeds_p84v2_v1` | `scripts/57_figS_comparator_seeds.py` | `pre_md_run/B_perseed_{stability,macro}.tsv`, `_partial/p84v2/31b_structcs_{45A,6A}/*.npz`, `validation_fix/F_structcs_{differences_20k,macro}.tsv`, `validation_fix/structcs_randomsplit_deterministic/*.csv` | Comparator subset (88 structures), all four cells: ensemble and the three per-seed values |
| `Figure_S7.png/.pdf` | `figS_legacy_p28_vs_p84v2_v1` | `scripts/57_figS_legacy_p28.py` | p28: `protocol_main.csv`, `benchmark_v2.csv`, `similarity_curve.csv`, `peptide_aware.csv`, `data/processed/pair_maps.npz` (shape check); p84v2: `stats_20k_seed42.tsv`, `percomplex_tcr70_p84v2.npz`, `percomplex_tcr70_paired_p84v2.tsv`, `benchmark_p84v2.csv`, `similarity_curve_p84v2.csv`, `peptide_aware_p84v2.csv` | Legacy 28-position headline results beside p84v2, with the outcome-informed selection disclosed |
| `Figure_S8.png/.pdf` | `figS_allele_trbv_p84v2_v1` | `scripts/57_figS_allele_trbv.py` | `percomplex_tcr70_paired_p84v2.tsv`, `pairwise_ours_p84v2.npz`, `protocol_main_p84v2.csv`, `complex_metadata.tsv` | Per-complex pairwise AUPRC under tcr70 by MHC allele and by TRBV (descriptive, no test) |

| `Figure_S9.png/.pdf` | `figS_sweep_composition_p84v2_v1` | `scripts/57_figS_sweep_composition.py`; composition by `scripts/57_sweep_composition_p84v2.py` | `similarity_curve_p84v2.csv`, `similarity_curve_composition_p84v2.csv`, `similarity_curve_composition_folds_p84v2.csv` (from `tcr_identity_matrix.npy` and the sweep's rng(0) folds) | Identity sweep (Figure 3B) with fold composition and the train/test similarity actually achieved |

**S4–S9 added 26 October 2026** (branch `si-figures`; `results/review/SI_FIGURES_REPORT.md`). S9's composition was
computed after the author authorised it; the `mean_max_sim` column of `similarity_curve*.csv` is an upper bound, not
the achieved similarity (see the report).

## Files in `figures/` that are NOT final figures

**Superseded: do NOT use** (p82v1 analysis, or pre-correction drafts):

| File | Why |
|---|---|
| `fig3_inflation_p82v1_v1.*`, `fig4_benchmark_p82v1_v1.*`, `fig5_footprint_p82v1_v2.*`, `fig6_ablation_p82v1_v1.*`, `figS_peptide_p82v1_v1.*`, `figS_benchmark_distributions_p82v1_v1.*` | p82v1 versions, superseded by p84v2 (residue representation corrected) |
| `fig3_inflation_v3.png`, `fig4_benchmark_v3.png`, `fig5_geometry_v4.png`, `fig6_peptide_failure_v3.png`, `fig2_dataset_final.png` | August drafts on the superseded 28-position candidate set |
| `figures/pymol/p82v1_*`, `figures/pymol/panel*.png`, `figures/pymol/session.pse`, `figures/pymol/cropped/*`, `figures/pymol/fig6_sites_p82v1.tsv` | p82v1 or August PyMOL intermediates |

**Classified as non-figure files:**

| File | What it is |
|---|---|
| `fig2_caption.txt` | The final Figure 2 caption as plain text (replaced 26 October 2026; the August text it held is gone). A copy: the source is `docs/manuscript/FIGURE_CAPTIONS.md` |
| `palette_reference.png/.pdf` | A swatch sheet of the `style.py` palette. A style reference, not a manuscript or SI figure |
| `figures/pymol/p84v2_*.png`, `p84v2_session.pse`, `fig6_sites_p84v2.tsv` | Current PyMOL panels and site tables, consumed by Figures 5 and 6. They are components, not figures |

**Not classifiable:** none. Every file in `figures/` is either a final figure, a component of one, a style reference,
or superseded.
