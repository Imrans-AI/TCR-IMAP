# Reproducing the reported results

Every command below runs from the **package root** with the environment of `env/environment.yml`
active. The first two groups need **no GPU, no network and no training**: they read stored arrays.
Retraining is a separate, optional section at the end.

```bash
conda env create -f env/environment.yml && conda activate tcr-imap
cd <package root>
```

If a command reports a missing array, you have the GitHub tree without the frozen data. See
`DATA.md`.

---

## 1. Start here: integrity, then the worked example

| command | what it checks or reproduces | expected output | cost |
|---|---|---|---|
| `python tests/test_frozen_artifacts.py` | cohort size, cluster counts, the five protocols, zero tcr70 cluster overlap, the indel-ratio similarity metric, label shape and sparsity, prediction/label keying, the approved Figure 1 hash, absence of absolute home paths | `8/8 checks passed` | seconds |
| `python scripts/worked_example.py` | the released arrays against the frozen labels, then the pairwise tcr70 comparisons and the identity sweep | five `[ok]` integrity lines, then the values below | < 1 min, CPU |

`worked_example.py` prints each computed value beside the value in the paper. It should reproduce:

```
TCR-IMAP vs pair-frequency prior   +0.0507 [+0.0215, +0.0841], p < 0.001, higher on 112/174
TCR-IMAP vs TEIM-Res               +0.0873 [+0.0592, +0.1149], p < 0.001, higher on 143/174
macro-AUPRC   TCR-IMAP 0.3598   prior 0.3090   TEIM-Res 0.2725
```

followed by the authoritative identity sweep (20,000 resamples, seed 42), whose 0.70 row is
**+0.0274, 95 % CI [+0.0011, +0.0542], p = 0.041**.

---

## 2. Tables and analyses from stored results

| command | reproduces | cost |
|---|---|---|
| `python scripts/59_verify_review_findings.py` | what each protocol actually holds out; the indel-distance check; the receptor-species composition | seconds |
| `python scripts/69_sweep_authoritative_stats.py` | the identity sweep at 20,000 resamples and seed 42, and the paired training-size effect | ~2 min |
| `python scripts/63_control_estimands.py` | the training-size-matched control's estimands and the paired contrast | ~1 min |
| `python scripts/64_control_diagnostics.py` | the control's similarity and training-size diagnostics | ~1 min |
| `python scripts/65_same22_comparator.py` | the comparator subset on a common set of 22 test complexes | seconds |
| `python scripts/66_fr4_sensitivity.py` | canonical-CDR-only rescoring across 14 comparisons | ~2 min |
| `python scripts/68_comparator_gate_sensitivity.py` | the comparator-fidelity gate sensitivity, gated and ungated | ~1 min |
| `python scripts/61_validation_split_paired.py` | the paired before/after difference for the validation-split correction | ~1 min |
| `python scripts/55_stats_p84v2.py` | the full statistics table behind Figures 4 and 6 and SI Table S1 | ~3 min |
| `python scripts/67_si_control_tables.py` | SI Tables S9 and S10 | seconds |
| `python docs/manuscript/si_tables/build_si_tables.py` | all SI tables, Markdown and CSV, from the committed tables | seconds |
| `python scripts/58_verify_captions.py` | every number in every figure caption against the tables the figure reads | seconds |

`58_verify_captions.py` must end with `ALL CAPTION NUMBERS VERIFIED`. It is the cheapest way to tell
whether anything in this package has drifted from what the captions assert.

---

## 3. Figures from stored results

All figure scripts write a 600-dpi PNG and a vector PDF next to each other under `figures/`, print a
panel-to-source validation table, and read only stored arrays.

| command | figure |
|---|---|
| `python scripts/38b_fig2_dataset_p82v1.py` | Figure 2, redundancy and leakage |
| `python scripts/55_fig3_inflation_p84v2.py` | Figure 3, identity sweep and docking-angle collapse |
| `python scripts/55_fig4_benchmark_p84v2.py` | Figure 4, benchmark against the comparators |
| `python scripts/55_fig5_footprint_p84v2.py` | Figure 5, observed contact footprint (PyMOL panels are included pre-rendered) |
| `python scripts/55_fig6_ablation_p84v2.py` | Figure 6, ablations and the two showcase complexes |
| `python scripts/56_fig7_md_occupancy.py` | Figure 7, MD occupancy of crystal-negative pairs |
| `python scripts/55_figS_peptide_p84v2.py` | Figure S1 |
| `python scripts/55_figS_benchmark_dist_p84v2.py` | Figure S2 |
| `python scripts/56_figS_md_characterisation.py` | Figure S3 |
| `python scripts/57_figS_validation_split.py` | Figure S4 |
| `python scripts/57_figS_cluster_macro.py` | Figure S5 |
| `python scripts/57_figS_comparator_seeds.py` | Figure S6 |
| `python scripts/57_figS_legacy_p28.py` | Figure S7 |
| `python scripts/57_figS_allele_trbv.py` | Figure S8 |
| `python scripts/57_figS_sweep_composition.py` | Figure S9 |

`figures/FIGURE_MANIFEST.md` maps each figure to its script and to the tables it reads.

**Figure 1 is an approved asset and is not regenerated.** Its editable source is
`figures/source/fig1_arial_20261102/design_figure.py`; `figures/main/Figure_1.SOURCE.md` records the
provenance, the revision history and the SHA-256 the manuscript builder enforces. Re-running the
source script reproduces the figure, but the raster is not byte-stable across matplotlib versions,
so a different matplotlib gives a visually identical figure with a different hash. The superseded
`scripts/70_fig1_workflow.py` refuses to run without `--allow-superseded` and can only write to
`figures/previous_versions/`.

Rebuilding the manuscript and the Supporting Information needs Node.js and `docx@9.15.0`:

```bash
NODE_PATH=<dir containing node_modules> node docs/manuscript/build_manuscript_docx_v3.js
NODE_PATH=<dir containing node_modules> node docs/manuscript/build_si_docx_v3.js
```

`docs/manuscript/BUILD_NOTES_v3.md` explains why the `docx` version is pinned: 9.7.1 silently breaks
all thirteen equations while leaving the XML text correct.

---

## 4. Optional: retraining

Nothing the paper concludes requires this. Retraining regenerates the stored predictions, which the
sections above consume. **A GPU is required.** All timings are measured on one
NVIDIA GeForce RTX 5070 Ti (16 GB), the machine every reported result was produced on; scale
accordingly.

| command | what it retrains | resources |
|---|---|---|
| `bash scripts/53_run_all_p84v2.sh` | the full p84v2 pipeline: pairwise and per-residue benchmarks, comparators, ablations, the sweep and the docking-angle regression | ≈ 10–14 h GPU, fits in 16 GB VRAM, ~55 GB disk for the resume shards |
| `python scripts/60_size_matched_control.py` | the training-size-matched control: 3 thresholds × 5 folds × 5 draws × 3 seeds = 225 fits | ≈ 10 min GPU |
| `python scripts/62_filtered_rerun_percomplex.py` | the per-complex arrays behind the authoritative sweep, 7 thresholds | ≈ 3 min GPU |
| `python scripts/19_train_pairmap.py` | one pairwise model, for inspection | minutes |

Both control scripts carry **validation gates** against the frozen aggregates and will refuse to
overwrite anything if a recomputed value disagrees. `scripts/53_run_all_p84v2.sh` takes a
single-instance lock so two runs cannot interleave.

Determinism: the stages re-run under the p84v2 corrections set `CUBLAS_WORKSPACE_CONFIG=:4096:8`
before importing torch, turn `cudnn.benchmark` off, set `cudnn.deterministic`, and call
`torch.use_deterministic_algorithms(True, warn_only=True)`. Even so, the 3 × 3 `Conv2d` weight
gradient in the TEIM-Res reimplementation is not bitwise reproducible across cuDNN builds; the
stored predictions are authoritative, and `results/review/TEIM_DETERMINISM_REPORT.md` documents it.

### Rebuilding the benchmark from structures

```bash
python scripts/01_fetch_stcrdab.py            # needs network; see data/SOURCES.md
python scripts/02_build_index.py
python scripts/22_build_metadata.py
python scripts/23_make_protocols.py           # WILL NOT overwrite the frozen protocols
python scripts/53_build_pair_labels_p84v2.py
python scripts/53_build_graphs_p84v2.py
```

The structures are third-party and are not redistributed; `data/SOURCES.md` and
`THIRD_PARTY_NOTICES.md` give the retrieval instructions. The frozen protocols and metadata are
shipped precisely so this step is not needed to reproduce anything reported.

### One retraining step needs a file that is not in the package

Re-running the comparator's **own single-split design** (`scripts/31b_structcs_setting_p82v1.py`,
which the p84v2 runner calls) reads `data/raw/sgtpmi/data/STRUCT-CS/PDBids.txt`, the test-set
identifier list from the released Struct-TCS repository. It is referenced rather than redistributed
because permission is unclear; `THIRD_PARTY_NOTICES.md` says where to get it. The **results** of
that design are included as `results/tables/structcs_setting_p84v2_{4p5A,6A}.csv`, so no reported
number depends on having the file.

### Molecular dynamics

Running the simulations needs GROMACS 2023.2 and the four CHARMM-GUI systems, which are not
redistributed; `md/systems/BUILD_MANIFEST.tsv` records how to rebuild them.
`md/MD_SYSTEM_SELECTION_CRITERIA.md` carries the predeclared selection, protocol and analysis
definitions. Twelve 200 ns replicates are roughly 1,500 GPU-hours in total. The occupancy analysis
itself (`md/analysis/`) runs on the trajectories; its outputs are already in `results/tables/md_*.csv`,
so the analysis can be checked without rerunning anything.
