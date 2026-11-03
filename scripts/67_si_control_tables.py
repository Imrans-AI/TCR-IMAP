#!/usr/bin/env python3
"""Supporting Information tables S9 and S10: the size-matched control's detailed memberships and
diagnostics, and the two sensitivities computed from stored artifacts.

Read-only: reads the tables written by scripts 60 and 62-66 and emits Markdown + CSV into
docs/manuscript/si_tables/, in the same style as build_si_tables.py. No training, no GPU.

usage: python scripts/67_si_control_tables.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/"docs/manuscript/si_tables"
RV = ROOT/"results/review"


def fmt(c, v):
    """Per-column formatting: thresholds and sizes are not 4-dp quantities."""
    if pd.isna(v): return ""
    if isinstance(v, (int,)) and not isinstance(v, bool): return str(v)
    if isinstance(v, float):
        if c == "threshold": return f"{v:g}"
        if c in ("training mean",): return f"{v:.1f}"
        if c in ("n", "clusters", "draw seed"): return f"{v:g}"
        return f"{v:.4f}".replace("-", "\u2212")
    return str(v).replace("-", "\u2212") if c == "setting" else str(v)


def md_table(df, floats=4):
    d = df.copy()
    for c in d.columns:
        d[c] = d[c].map(lambda v, c=c: fmt(c, v))
    head = "| " + " | ".join(d.columns) + " |"
    sep = "|" + "|".join(["---"]*len(d.columns)) + "|"
    body = "\n".join("| " + " | ".join(r) + " |" for r in d.values.tolist())
    return "\n".join([head, sep, body])


# ------------------------------------------------------------------ S9
est = pd.read_csv(RV/"control_estimands/control_estimands.tsv", sep="\t")
diag = pd.read_csv(RV/"control_estimands/similarity_diagnostics.tsv", sep="\t")
mem = pd.read_csv(RV/"size_matched_control/memberships.tsv", sep="\t")
drw = pd.read_csv(RV/"size_matched_control/per_draw_metrics.tsv", sep="\t")

e = est[["threshold", "principal_mean_of_draw_differences", "principal_lo", "principal_hi",
         "draw_min", "draw_max", "draw_range", "filtered_delta", "filtered_lo", "filtered_hi",
         "contrast_filtered_minus_sizematched", "contrast_lo", "contrast_hi",
         "secondary_ensemble_delta", "secondary_ensemble_lo", "secondary_ensemble_hi"]].copy()
e.columns = ["threshold", "principal Δ_rand", "lo", "hi", "draw min", "draw max", "draw range",
             "Δ_sim", "lo ", "hi ", "paired contrast", "lo  ", "hi  ",
             "ensemble Δ_rand (secondary)", "lo   ", "hi   "]

d9 = diag[["regime", "threshold", "draw_seed", "train_sizes_per_fold", "train_size_mean",
           "max_over_test", "mean_of_fold_maxima", "mean_per_test_nn"]].copy()
d9.columns = ["regime", "threshold", "draw seed", "training per fold", "training mean",
              "max over all test complexes", "mean of fold maxima", "mean per-test nearest-training"]

pd_tab = drw.copy()
pd_tab.columns = ["threshold", "draw seed", "n", "model", "prior", "Δ", "lo", "hi", "p"]

s9 = f"""# Table S9. Training-size-matched control: estimands, diagnostics and memberships

**Table S9.** The control specified in `results/review/SIZE_MATCHED_CONTROL_PLAN.md`, frozen before any model was
trained, and reported in `SIZE_MATCHED_CONTROL_RESULT.md`. For each threshold and fold the number of training
complexes retained by the similarity filter was drawn instead at random from the unrestricted pool, with no
similarity filter, on the same test complexes; five draws with seeds fixed in advance (1001–1005), three training
seeds each, 225 fits in total.

## S9a. Estimands

The **principal** estimate is the mean over draws of the per-complex difference between the model and the pair-frequency prior, bootstrapped over the 69 clu70 clusters (20,000 resamples, seed 42), and is **conditional on the five realized training draws**.
The **draw range** is a five-point spread over a shared test set and is **not a confidence interval**.
The **paired contrast** is computed per complex on the same 174 test complexes under the same cluster resamples.
The **ensemble** quantity averages predictions across draws before scoring; because average precision is nonlinear it is a different estimand and is reported as **secondary**.
All p values fall below the 1/20,000 = 5 × 10⁻⁵ floor and are reported as p < 0.001.

{md_table(e)}

## S9b. Similarity and training-size diagnostics

Three distinct quantities, kept separate and never conflated.
The **maximum over all evaluated test complexes** of a test complex's similarity to its fold's training set.
The **mean of the five fold maxima**, which is smaller.
The **mean over test complexes of each complex's nearest-training similarity**, smaller again.
Training sizes are exact per fold; any mean is labelled as a mean.

{md_table(d9)}

## S9c. Per-draw results

{md_table(pd_tab)}

Full subset memberships — the training PDB identifiers drawn for every (threshold, fold, draw) — are in
`results/review/size_matched_control/memberships.tsv` ({len(mem)} subsets) and in the deposit.
"""

# ------------------------------------------------------------------ S10
s22 = pd.read_csv(RV/"same22/same22_comparison.tsv", sep="\t")
t10a = s22[["cutoff", "comparator", "setting", "n", "n_clusters", "ours", "comp", "diff", "lo", "hi", "p"]].copy()
t10a.columns = ["cutoff", "comparator", "setting", "n", "clusters", "TCR-IMAP", "comparator AUPRC",
                "Δ", "lo", "hi", "p"]

# 66_fr4_sensitivity.py was extended from the three sweep thresholds to every comparison whose
# predictions are stored, and its schema changed with it: the columns are now
# group/comparison/scope/n/a/b/diff/lo/hi/p. Read what the file actually has.
fr4 = pd.read_csv(RV/"fr4_sensitivity/fr4_sensitivity.tsv", sep="\t")
assert list(fr4.columns) == ["group", "comparison", "scope", "n", "a", "b",
                             "diff", "lo", "hi", "p"], list(fr4.columns)
t10b = fr4[["group", "comparison", "scope", "n", "a", "b", "diff", "lo", "hi", "p"]].copy()
t10b.columns = ["evaluation", "comparison", "scope", "n", "TCR-IMAP", "comparator",
                "Δ", "lo", "hi", "p"]
chg = pd.read_csv(RV/"fr4_sensitivity/fr4_change_summary.tsv", sep="\t")
n_comparisons = len(chg)
worst = chg.iloc[chg["change"].abs().idxmax()]
# U+2212 for the sign, as everywhere else in the figures and tables.
worst_change = f"{worst['change']:+.4f}".replace("-", "−")
teim = pd.read_csv(RV/"fr4_sensitivity/fr4_teim.tsv", sep="\t")

s10 = f"""# Table S10. Common-test-set and target-definition sensitivities

**Table S10.** Two sensitivities computed entirely from stored prediction and per-complex arrays. Nothing was
retrained for either.

## S10a. The comparator subset on a common test set

The reported comparator-subset contrast sets a single random split of 22 test complexes against tcr70 with 88, so it
confounds the partitioning protocol with the set of complexes scored. The 22 are a subset of the 88, so the tcr70
evaluation can be restricted to exactly those 22. The 88-complex rows are preserved and remain the reported result.
Paired cluster bootstrap over clu70, 20,000 resamples, seed 42.

{md_table(t10a)}

At 6 Å the sign of the difference from Struct-TCS **reverses between protocols on identical test complexes** (−0.0949 under the single random split against +0.0735 under tcr70).
That locates the instability in the partitioning design rather than in the choice of test complexes.
Seed-to-seed variation within a fixed split (Figure S6) does not measure variation between splits and is not used to infer it.

## S10b. Canonical-CDR-only scoring, FR4 excluded

The 75-row TCR representation keeps 14 J-region framework rows (IMGT 118–124 of both chains).
Those rows carry **0 of the 5,281 contacts** at 4.5 Å but **204,624 of the 893,928 resolved pairs** — 22.9 % of the evaluated negatives and none of the positives.
Rescoring the stored predictions over canonical CDR rows only (IMGT ≤ 117) tests how much each comparison depends on them.
Paired cluster bootstrap over clu70, 20,000 resamples, seed 42. The *TCR-IMAP* and *comparator* columns are macro-AUPRC.

{md_table(t10b)}

**What this covers.** {n_comparisons} comparisons could be rescored: TCR-IMAP against the pair-frequency prior on the
pairwise task under all five protocols, TCR-IMAP against TEIM-Res under tcr70, TCR-IMAP against Struct-TCS and against
Seq-TCS on the per-residue task over the 88-complex comparator subset, and the three sweep thresholds in both the
similarity-filtered and the size-matched-random regime. Across those {n_comparisons} the largest change in the paired
difference is **{worst_change}** ({worst['group']}, {worst['comparison']}); every other comparison moves by at most
0.0002, and no interval crosses zero that did not already. TEIM-Res’s own macro-AUPRC under tcr70 moves from
{teim.macro_auprc.iloc[0]:.4f} to {teim.macro_auprc.iloc[1]:.4f}. Within this listed set the FR4 rows change no verdict,
which is a statement about these {n_comparisons} comparisons and not about every result in the paper.

**What this does not cover.** The **predeclared primary comparison** — per-residue macro-AUPRC under tcr70, TCR-IMAP
versus Struct-TCS, over the full 174-complex benchmark — **cannot be rescored at all**. That runner stores per-complex
metrics and no prediction array, and writes no checkpoint, so there is nothing to restrict to canonical CDR rows; the
files under `results/checkpoints/` are from the superseded 28-position candidate set under partitions that are not the
frozen protocols, and are not the reported runs (`results/review/FR4_CHECKPOINT_AUDIT.md`). Completing this one
comparison would require retraining the per-residue benchmark with predictions persisted. It is reported as
outstanding, **not** estimated from the rows above: the 88-complex comparator-subset figures are a different
evaluation on a different set of complexes and are not a proxy for it. The gap is in one robustness check; the primary
comparison itself was computed over its full predeclared candidate universe.
"""

(OUT/"S9_size_matched_control.md").write_text(s9.rstrip() + "\n")
(OUT/"S10_sensitivities.md").write_text(s10.rstrip() + "\n")
est.to_csv(OUT/"S9_control_estimands.csv", index=False)
diag.to_csv(OUT/"S9_control_diagnostics.csv", index=False)
drw.to_csv(OUT/"S9_control_per_draw.csv", index=False)
mem.to_csv(OUT/"S9_control_memberships.csv", index=False)
s22.to_csv(OUT/"S10_same22.csv", index=False)
fr4.to_csv(OUT/"S10_fr4_sensitivity.csv", index=False)
print("wrote S9_size_matched_control.md, S10_sensitivities.md and 6 CSVs")
