# Table S5. Identity-threshold sweep: composition and results

**Table S5.** For each threshold on the TCR Levenshtein similarity ratio, the training complexes kept in each of the
five folds, the train/test similarity actually achieved, and pairwise macro-AUPRC of TCR-IMAP and the pair-frequency
prior with their paired difference.

| Threshold | n test | n training per fold (folds 0–4) | Mean n training | Achieved max similarity | Achieved mean similarity | `mean_max_sim` column (upper bound) | TCR-IMAP | Prior | Δ | 95 % CI | p | 95 % CI (superseded) | p (superseded) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| none | 174 | 139 / 139 / 139 / 139 / 140 | 139.2 | 1.000 | 0.890 | 0.8977 | 0.6313 | 0.3563 | +0.2750 | [+0.2144, +0.3203] | < 0.001 | [+0.2155, +0.3202] | < 0.001 |
| < 0.95 | 174 | 91 / 94 / 96 / 95 / 94 | 94.0 | 0.948 | 0.719 | 0.7482 | 0.3949 | 0.3172 | +0.0777 | [+0.0336, +0.1166] | < 0.001 | [+0.0347, +0.1172] | < 0.001 |
| < 0.90 | 174 | 88 / 90 / 93 / 94 / 90 | 91.0 | 0.862 | 0.699 | 0.7189 | 0.3675 | 0.3155 | +0.0520 | [+0.0177, +0.0841] | 0.003 | [+0.0194, +0.0839] | 0.005 |
| < 0.85 | 174 | 86 / 86 / 88 / 93 / 86 | 87.8 | 0.841 | 0.682 | 0.7028 | 0.3526 | 0.3098 | +0.0428 | [+0.0101, +0.0747] | 0.009 | [+0.0100, +0.0755] | 0.009 |
| < 0.80 | 174 | 84 / 80 / 86 / 92 / 84 | 85.2 | 0.790 | 0.672 | 0.6889 | 0.3406 | 0.3066 | +0.0340 | [+0.0063, +0.0609] | 0.013 | [+0.0053, +0.0609] | 0.016 |
| < 0.75 | 174 | 84 / 80 / 83 / 92 / 84 | 84.6 | 0.740 | 0.671 | 0.6861 | 0.3393 | 0.3063 | +0.0330 | [+0.0044, +0.0605] | 0.021 | [+0.0034, +0.0607] | 0.028 |
| < 0.70 | 174 | 72 / 69 / 76 / 78 / 73 | 73.6 | 0.699 | 0.649 | 0.6627 | 0.3271 | 0.2997 | +0.0274 | [+0.0011, +0.0542] | 0.041 | [−0.0001, +0.0543] | 0.052 |

**Notes.**
- The sweep uses one fixed random five-fold partition (`np.random.default_rng(0)`, independent of the frozen
  protocols; folds of 35, 35, 35, 35 and 34). At each threshold every non-test complex whose similarity to **any** test
  complex of its fold reaches the threshold is removed from that fold's training set; a fold is used only if at least 30
  training complexes remain (all five are used at every threshold). The test composition never changes: all 174
  complexes are evaluated at every threshold. The model and prior are retrained at every threshold.
- Achieved similarity: for each test complex, its maximum similarity to the training complexes its fold kept; the
  maximum and the mean over the 174 test complexes.
- **The `mean_max_sim` column recorded with the sweep is not the achieved similarity.** For each test complex it takes
  the maximum similarity to any *other complex of the 174* whose similarity is below the threshold, including complexes
  in the same test fold and training complexes the filter removed, and averages over test complexes. It is therefore
  an upper bound on the achieved mean (here 0.008–0.029 above it, and different even with no restriction, because it
  also counts same-fold neighbors).
- Δ = TCR-IMAP − prior. The reported **95 % CI and p use the project standard: paired clu70 cluster bootstrap,
  20,000 resamples, seed 42**, p floored at 1/20,000 and reported as < 0.001 at the floor.
- The two rightmost columns carry the **superseded** values first recorded for this stage, which used the
  2,000-resample seed-0 default. The sweep originally stored no per-complex arrays, so it could not be recomputed at
  the standard; those arrays were later regenerated under the stage's own configuration and seeds, each threshold
  reproducing the frozen aggregate to four decimal places, which is what allows the recomputation.
  **The point estimates are identical under both settings** (asserted in the build). The intervals and p values differ
  only by Monte Carlo error in the bootstrap, which 2,000 resamples do not resolve near an interval endpoint; at
  < 0.70 that moves the interval endpoint across zero (p 0.052 → 0.041). This is a numerical property of the
  resampling, not a second analysis and not evidence about evaluation design.
- Similarity uses the clustering strings (see Table S4).

**Sources.** `results/tables/similarity_curve_p84v2.csv` (results, `mean_max_sim`),
`results/tables/similarity_curve_composition_p84v2.csv` and `similarity_curve_composition_folds_p84v2.csv` (composition;
built by `scripts/57_sweep_composition_p84v2.py` from `tcr_identity_matrix.npy` and the sweep's own fold seed;
`results/review/SI_FIGURES_REPORT.md`, section S9). Sweep procedure: `scripts/25b_similarity_curve_p82v1.py`, executed on
the p84v2 inputs by `scripts/53_similarity_curve_p84v2.py`.
