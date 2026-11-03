# Table S2. Ablations of TCR-IMAP (tcr70, pairwise task)

**Table S2.** Change in pairwise macro-AUPRC under tcr70 when one input or architectural component is removed, for all
eight ablations. Δ = full model − ablated model (positive, the ablation lowered performance). Full TCR-IMAP:
macro-AUPRC 0.3598, n = 174 complexes, 69 clusters.

| Ablation | Code | Ablated macro-AUPRC | Δ | 95 % CI | 99.375 % simultaneous interval | p (raw) | p (Holm, 8) | Survives Holm |
|---|---|---|---|---|---|---|---|---|
| MHC residue identity masked | `no_mhc` | 0.3418 | +0.0180 | [+0.0008, +0.0346] | [−0.0054, +0.0408] | 0.040 | 0.160 | no |
| TCR residue identity masked | `no_tcr` | 0.3261 | +0.0337 | [+0.0044, +0.0662] | [−0.0052, +0.0792] | 0.021 | 0.104 | no |
| TCR residues shuffled within each loop | `shuffle_tcr` | 0.3095 | +0.0503 | [+0.0265, +0.0759] | [+0.0189, +0.0859] | < 0.001 | < 0.001 | yes |
| MHC residues shuffled across positions | `shuffle_mhc` | 0.3314 | +0.0283 | [+0.0103, +0.0458] | [+0.0040, +0.0524] | 0.002 | 0.014 | yes |
| CDR3 + FR4 tokens removed | `no_cdr3` | 0.3477 | +0.0121 | [−0.0044, +0.0293] | [−0.0103, +0.0369] | 0.161 | 0.483 | no |
| Germline CDR1/CDR2 tokens removed | `no_germline` | 0.3137 | +0.0461 | [+0.0293, +0.0639] | [+0.0228, +0.0708] | < 0.001 | < 0.001 | yes |
| Cross-attention removed | `no_cross` | 0.3606 | −0.0008 | [−0.0111, +0.0104] | [−0.0148, +0.0146] | 0.875 | 1.000 | no |
| Bilinear head replaced by dot product | `dot_head` | 0.3636 | −0.0039 | [−0.0170, +0.0104] | [−0.0214, +0.0169] | 0.581 | 1.000 | no |

**Notes.**
- Each ablated model was trained and evaluated exactly as the full model (same folds, epochs, optimizer and seeds 42–44,
  predictions averaged), on the same 174 complexes.
- 95 % CI, simultaneous interval and raw p come from the same 20,000 cluster-bootstrap resamples (pooled clu70, seed 42),
  **recomputed by `build_si_tables.py` from the stored per-complex arrays** and asserted equal to the committed values.
  The simultaneous interval is the Bonferroni 1 − 0.05/8 = 99.375 % percentile interval (quantiles 0.3125 % and 99.6875 %).
- Holm step-down adjustment across the eight ablations. "Survives Holm": Holm-adjusted p < 0.05; in every row this
  coincides with the simultaneous interval excluding zero (asserted).
- p is floored at 1/20,000 and reported as < 0.001 at or below 0.0005. The `p` column of `ablations_p84v2.csv` is the
  run's own 2,000-resample value.
- CDR3 + FR4: the CDR3 token block spans IMGT 105–124, including the J-region framework positions 118–124.

**Sources.** `results/tables/_partial/p84v2/29b_ablations/{full,no_mhc,no_tcr,shuffle_tcr,shuffle_mhc,no_cdr3,no_germline,no_cross,dot_head}.npz` (per-complex AUPRC); clusters
from `results/tables/pairwise_ours_p84v2.npz`; means checked against `results/tables/ablations_p84v2.csv`; all
statistics checked against `results/review/figures_p84v2/stats_20k_seed42.tsv`. Ablation definitions:
`scripts/29b_ablations_p82v1.py` (`apply_ablation`, `NoCross`, `DotHead`), executed on the p84v2 inputs by
`scripts/53_ablations_p84v2.py`.
