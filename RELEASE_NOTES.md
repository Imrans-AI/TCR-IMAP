# Release notes — 1.0.0-rc1

**Status: release candidate, prepared locally for author review. Not published, not pushed, not
deposited.** No repository URL and no DOI exist. See `PUBLICATION_CHECKLIST.md`.

## Versioning

This is the first version, so there is no prior scheme to preserve. The proposed scheme is
**semantic versioning**, which is what Zenodo's `version` field expects:

* `1.0.0-rc1` — this candidate: everything assembled and tested locally, awaiting author review.
* `1.0.0` — the version deposited and tagged when the manuscript is accepted, with the licence
  decided and the repository and archive records created.
* `1.0.x` — documentation and packaging fixes that change no result.
* `1.1.0` — new analyses or data that do not invalidate the frozen protocols or the stored
  predictions.
* `2.0.0` — any change to the frozen protocols, the candidate set or the contact definition. These
  define the benchmark; changing one of them is a new benchmark, not a patch.

**No tag has been created.** `PACKAGE_INFO.json` records the source revision the packages were
staged from, which is what identifies this candidate until a tag exists.

## What is in this candidate

The frozen apparatus behind the reported results: 174 curated TCR–peptide–MHC class I complexes with
their contact labels and candidate masks, five partitioning protocols fixed before the results were
produced, the complete stored pairwise predictions for all five protocols, per-complex per-residue
metrics, the within-fold frequency priors, the comparator reimplementations with a deviation table,
the training-size-matched control with the plan that was committed before any downsampled model was
trained, the authoritative identity sweep, the canonical-CDR sensitivity, the molecular-dynamics
analysis code and its twelve output tables, every final figure with its source, and a worked example
that reproduces the headline pairwise numbers on a CPU in under a minute.

## Notable corrections carried in this candidate

| what | correction |
|---|---|
| identity sweep uncertainty | the intervals and p values are the authoritative **20,000 resamples at seed 42**, not the function defaults. The earlier `[−0.0001, +0.0543], p = 0.052` differed from `[+0.0011, +0.0542], p = 0.041` by **Monte Carlo error alone** from an identical point estimate of +0.0274. Both are shipped, labelled, in SI Table S5 |
| the `double` protocol | identified as an **intersection-selected legacy partition, not a double holdout**: its held-out groups remain in training. Relabelled "intersection (legacy)", removed from any strictness ordering; the frozen key is unchanged |
| the allele protocol | described as **holdout by available allele annotation**: 108 of 174 complexes carry no allele annotation, and 32 of 64 evaluated complexes have an identical 84-residue MHC window sequence in training |
| receptor species | "174 human TCR complexes" replaced by the measured composition: human 165, mouse 8, mixed chimpanzee/human α 1. All 174 have a human MHC class I heavy chain |
| the similarity metric | `Levenshtein.ratio` normalises the **indel** distance, so a substitution costs 2. The definition, the operation costs and the software citation are corrected, and a test asserts it |
| comparator fidelity | a full provenance and deviation table per comparator, and the statement that this is an architecture comparison under common inputs, **not** a reproduction of published performance |
| Figure 1 | the author-approved two-panel figure, revised the same day to Arial typography to match Figures 2–7. The builder verifies its SHA-256 and aborts on a mismatch |

## Known limitations shipped with this candidate

1. **The predeclared primary comparison is unresolved** at this sample size: per-residue macro-AUPRC
   under tcr70, TCR-IMAP versus Struct-TCS, +0.0312 [−0.0049, +0.0650], p = 0.088. Reported as it
   stands; no secondary comparison was promoted in its place.
2. **The canonical-CDR-only sensitivity does not cover that primary comparison.** It covers 14
   comparisons, largest change 0.0014. The per-residue benchmark writes **no checkpoint at all** and
   stores metrics rather than predictions, so the primary cannot be rescored without retraining.
   `results/review/FR4_CHECKPOINT_AUDIT.md` carries the evidence: `torch.save` appears zero times in
   the per-residue path, and the files under `results/checkpoints/` are from the superseded
   28-position candidate set under partitions that are not the frozen protocols. **No historical
   per-residue checkpoint or prediction array exists, and this package does not claim otherwise.**
3. **No sequence-to-input builder.** Applying the model outside the cohort needs an IMGT numbering
   step that is documented here but not implemented.
4. **Trajectories are not included** and no access arrangement is promised.
5. **No licence is assigned yet** — `LICENSING_DECISION_REQUIRED.md`.
