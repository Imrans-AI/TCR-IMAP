# Data dictionary

Every released array and table, field by field. Shapes are given for the 174-complex cohort with 75
TCR rows and 84 MHC columns. `FILE_MANIFEST.tsv` lists the size and SHA-256 of each file;
`DATA.md` says which of them are in the GitHub tree and which are in the Zenodo archive.

Three conventions hold throughout:

* **Complex order is fixed.** Every array indexed by complex uses the same order as
  `data/processed/complex_metadata.tsv`, so `pdb[i]` is the same complex in every file. The checks
  in `tests/test_frozen_artifacts.py` assert this.
* **Residues are keyed by (IMGT number, insertion code)**, never by IMGT number alone, and tokens
  are placed at the row or column their position assigns rather than packed from the start of a
  segment.
* **A pair is scored only where it is resolved.** `S = 1` marks pairs in which both residues are
  resolved in the crystal structure; `S = 0` cells are masked in training and in scoring and carry
  no meaningful prediction.

---

## Benchmark definition

### `data/processed/complex_metadata.tsv`: the cohort, 174 rows

| column | type | meaning |
|---|---|---|
| `pdb` | 4-char | PDB identifier, lowercase; the join key for every other file |
| `allele` | str | MHC class I allele as annotated by STCRDab; **empty for 108 complexes**, which is why the allele protocol is "holdout by available allele annotation" |
| `trav`, `trbv` | str | TCR α and β V-gene families |
| `peptide` | str | peptide description from the STCRDab summary |
| `date`, `dep_date` | date | PDB deposition date; `dep_date` defines the temporal protocol |
| `resolution` | float | crystallographic resolution, Å |
| `study_id` | str | author list of the originating study, used to detect same-study families |
| `docking_angle` | float | TCR docking angle in degrees, the AngleNet regression target |
| `cdr1a`…`cdr3b` | str | the six CDR loop sequences, one residue per occupied IMGT position |
| `tcr_seq` | str | the concatenated CDR slot sequence used for similarity and clustering |
| `cdr3_pair` | str | `CDR3α|CDR3β`, the key for duplicate families |
| `dup_family` | int | duplicate-family identifier: identical CDR3 pair and peptide |
| `clu70`, `clu80`, `clu90` | int | single-linkage TCR clusters at similarity ratio 0.70, 0.80 and 0.90, **69, 82 and 89 clusters**. `clu70` is the bootstrap resampling unit |

The similarity ratio is `Levenshtein.ratio` over `tcr_seq`, which normalises the **indel** distance:
a substitution costs 2, not 1. `tests/test_frozen_artifacts.py` asserts this, because a package that
normalised the Levenshtein distance instead would change every cluster and every split.

### `data/splits_v2/protocols.json`: the five frozen partitions

A mapping from protocol name to a list of folds, each fold a list of PDB identifiers. Folds with
fewer than three test complexes are dropped at evaluation, which is why the number of complexes
evaluated differs from 174.

| protocol | folds (test sizes) | evaluated | what it withholds |
|---|---|---|---|
| `random` | 5 (35, 35, 35, 35, 34) | 174 | nothing, the reference point |
| `tcr70` | 5 (39, 33, 25, 32, 45) | 174 | whole TCR clusters at ratio 0.70; **the primary protocol** |
| `allele` | 5 (26, 20, 15, 2, 3) | 64 | MHC allele groups, by available annotation; the 2-complex fold is below the minimum |
| `temporal` | 1 (60) | 60 | complexes deposited after a cut-off date |
| `double` | 3 (16, 5, 9) | 30 | **legacy.** Selects the intersection of held-out cluster and allele groups but trains on all non-test complexes, so its held-out groups remain in training. Displayed as "intersection (legacy)"; the key stays `double` so frozen files are not renamed |

### `data/processed/pair_maps_p84v2.npz`: labels, masks and the representation

| key | shape | dtype | meaning |
|---|---|---|---|
| `Y` | 174 × 75 × 84 | int8 | 1 where the residue pair is a contact: minimum inter-residue heavy-atom distance ≤ 4.5 Å. **5,281 contacts** |
| `S` | 174 × 75 × 84 | int8 | 1 where both residues are resolved. **893,928 resolved pairs**, so the contact rate is 0.59 % |
| `pdb` | 174 | str | complex order |
| `row_label` | 75 | str | TCR row label, e.g. `cdr3a_111A` |
| `row_segment` | 75 | str | one of `cdr1a cdr2a cdr3a cdr1b cdr2b cdr3b` |
| `row_imgt_number` | 75 | int | IMGT position |
| `row_insertion_code` | 75 | str | insertion code, empty for base positions; **7 CDR3 insertion positions occupy rows of their own** |
| `row_imgt_region` | 75 | str | `CDR1` (13 rows), `CDR2` (15), `CDR3` (33), `FR4` (14) |
| `col_label`, `col_imgt_number`, `col_insertion_code` | 84 | n/a | the MHC window: IMGT 50–90 of α1 (**41 columns**) and 1050–1090 of α2 (**43**), including the insertion columns 1061A and 1072A |
| `tcr_aa`, `mhc_aa` | 174 × 75, 174 × 84 | str | the amino-acid identity in each row/column; `-` where unoccupied. **These are the model inputs** |
| `tcr_resid`, `mhc_resid` | 174 × 75, 174 × 84 | str | the source residue identifier in the structure |
| `cutoff` | scalar | float | 4.5 |
| `candidate_set` | scalar | str | `p84v2` |
| `provenance` | scalar | str | free text recording how the file was built |

### `data/processed/tcr_graphs_p84v2.npz`: residue graphs for Struct-TCS

| key | shape | meaning |
|---|---|---|
| `A` | 174 × 75 × 75 | TCR residue adjacency; the Cβ/Cα construction described in `results/review/COMPARATOR_PROVENANCE.md`, which also records how it differs from the published atom definitions |
| `pdb`, `row_label` | n/a | the same keying as the label file |

### `data/interim/stcrdab_index.tsv`

The STCRDab summary rows for the cohort, including `alpha_organism` and `beta_organism`, the fields
the receptor-species composition is derived from (human in 165, mouse in 8, mixed chimpanzee/human α
in 1). `data/interim/vgene_cdr12.json` maps V genes to germline CDR1/CDR2 sequences;
`data/interim/prior_complexes.tsv` lists the complexes the descriptive prior is fitted on.

---

## Stored predictions and metrics

### `results/tables/pairmaps_{protocol}_p84v2.npz`: pairwise predictions, one per protocol

| key | shape | meaning |
|---|---|---|
| `pdb`, `clu70` | 174 | complex order and the bootstrap unit |
| `covered` | 174 | bool: did this protocol evaluate this complex. `random` 174, `tcr70` 174, `allele` 64, `temporal` 60, `double` 30 |
| `Y`, `S` | 174 × 75 × 84 | copies of the frozen labels and masks, so each file is self-checking |
| `pred_model` | 174 × 75 × 84 | TCR-IMAP's out-of-fold probability for every candidate pair |
| `pred_freq` | 174 × 75 × 84 | the within-fold pair-frequency prior, fitted on training folds only |

**Rows outside `covered` are NaN by design.** Range and finiteness only mean anything inside the
coverage mask, which is what the tests assert.

### `results/tables/percomplex_{protocol}_p84v2.npz`: per-residue metrics

`pdb`, `clu70`, `allele`, then `{method}_auprc` and `{method}_auc` for `Ours`, `Struct-TCS`,
`Seq-TCS`, `PosFreq` and `NN`, each 174 long. These are **metrics, not predictions**: the
per-residue benchmark stores aggregates per complex and no prediction array, which is exactly why
the canonical-CDR-only sensitivity cannot be computed for the predeclared per-residue primary
(`results/review/FR4_CHECKPOINT_AUDIT.md`). **No per-residue checkpoint or prediction array exists.**

### Other prediction and metric files

| file | contents |
|---|---|
| `results/tables/pairwise_ours_p84v2.npz` | per-complex pairwise AUPRC and AUC for TCR-IMAP (metrics only) |
| `results/tables/teim_res_pairwise_p84v2.npz` | TEIM-Res per-complex metrics **and** its full 174 × 75 × 84 prediction array |
| `results/tables/_partial/p84v2/31b_perseed/*.npz` | comparator-subset per-residue predictions, per training seed, at the 4.5 Å and 6 Å contact definitions |
| `results/tables/_partial/p84v2/25b_similarity_curve/thr_*.npz` | per-threshold arrays behind the identity sweep |
| `results/tables/_partial/p84v2/29b_ablations/*.npz` | per-ablation arrays behind Figure 6 |

### Aggregates

| file | contents |
|---|---|
| `benchmark_p84v2.csv`, `protocol_main_p84v2.csv` | the headline benchmark by method and protocol |
| `similarity_curve_p84v2.csv` | the sweep as originally recorded (2,000 resamples, seed 0) |
| `similarity_curve_p84v2_authoritative.csv` | **the authoritative sweep**: 20,000 resamples, seed 42. Point estimates are identical to the file above; the intervals and p values differ by Monte Carlo error alone, and only the 0.70 row changes verdict |
| `similarity_curve_composition_p84v2.csv`, `..._folds_p84v2.csv` | how many training complexes each threshold retained, per fold |
| `ablations_p84v2.csv` | ablation results with Holm-adjusted intervals |
| `cluster_macro_p84v2.csv` | the cluster-weighted sensitivity |
| `percomplex_maxsim.tsv` | each test complex's maximum similarity to any training complex, per protocol |
| `leakage_audit_dual.tsv` | duplicate families split between train and test, per protocol and fold |
| `tcr_identity_matrix.npy` | 174 × 174 similarity ratios over `tcr_seq` |
| `angle_runs_p84v2.jsonl` | one record per docking-angle run |
| `contact_freq_gdomain182_descriptive_p84v2.npz`, `gdomain182_region_classification_p84v2.csv`, `germline_diagonal_p84v2.csv` | the model-free contact footprint behind Figure 5 |
| `structcs_setting_p84v2_{4p5A,6A}.csv` | the comparator's own single-split design, at both contact definitions |
| `peptide_aware_p84v2.csv` | the peptide-input ablation behind Figure S1 |
| `md_*.csv` | the twelve molecular-dynamics output tables: occupancy by stratum, top-k, the docked-frame check, ρ values, RMSF, hydrogen bonds, characterisation and the run manifests |

### `results/review/`: the analyses, with their inputs and dated plans

| path | contents |
|---|---|
| `SIZE_MATCHED_CONTROL_PLAN.md` / `_RESULT.md` | the control's plan, **committed before any downsampled model was trained**, and its outcome |
| `size_matched_control/` | memberships, achieved similarities, per-draw and per-complex metrics, and the prediction arrays for the three thresholds |
| `control_estimands/` | the estimands and the paired contrast |
| `sweep_authoritative/` | the authoritative sweep statistics and the paired training-size effect |
| `filtered_rerun/` | per-complex arrays and predictions for the seven sweep thresholds, with a validation gate against the frozen CSV |
| `SWEEP_UNCERTAINTY_RECONCILIATION.md` | why two interval estimates differed: Monte Carlo error in the bootstrap, from an identical point estimate |
| `COMPARATOR_PROVENANCE.md` | published versus implemented, per comparator: features, graph construction, architecture, training, the execution fixes and the unavailable components |
| `COMPARATOR_FIDELITY_SENSITIVITY.md`, `comparator_fidelity/` | the gated and ungated variants and the paired contrast between them |
| `fr4_sensitivity/`, `FR4_CHECKPOINT_AUDIT.md` | the canonical-CDR-only rescoring over 14 comparisons, and the evidence that the predeclared per-residue primary cannot be rescored without retraining |
| `same22/` | the comparator subset on a common set of 22 test complexes |
| `V3_PART_A_VERIFICATION.md`, `V3_PART_C_VALIDATION_SPLIT.md`, `V3_FINAL_REPORT.md` | the external-review verification, the validation-split interpretation and the revision record |
| `CORRECTION_PLAN_P84V2.md` | the dated correction plan that restated the primary metric, protocol and comparison |
| `figures_p84v2/stats_20k_seed42.tsv` | the 20,000-resample statistics the figure captions are checked against |

---

## Figures

`figures/main/Figure_{1..7}.{png,pdf}` and `figures/supplementary/Figure_S{1..9}.{png,pdf}`, all
7.2 in wide with Arial embedded and no Type 3 font. `figures/FIGURE_MANIFEST.md` maps each figure to
its generating script and to the tables it reads. `figures/pymol/` holds the pre-rendered structural
panels used by Figures 5 and 6. `figures/source/fig1_arial_20261102/` is the editable source of the
approved Figure 1, and `figures/main/Figure_1.SOURCE.md` its provenance and revision history.
