# Audit record: what was included, what was left out, and why

Both packages are staged from the **current working tree** by
`release/packaging/stage_release.py`, not copied from any earlier snapshot. The staging rules are in
that script, in one declarative table, so what is in the package can be read off the code rather
than inferred from the result. `PACKAGE_INFO.json` records the source revision;
`FILE_MANIFEST.tsv` and `CHECKSUMS.sha256` record every file with its size and SHA-256.

The working repository is 14 GB. Publishing it indiscriminately would ship 917 MB of third-party
structures, 11 GB of molecular-dynamics trajectories and build bundles, 160 MB of superseded
checkpoints and a decade of intermediate state. The rules below are the filter.

## How `scripts/` was chosen

Not by copying the directory: of its 135 files, most belong to the superseded p28 and p82v1
generations. The package instead contains a **dependency closure**, computed at staging time, over
every command the documentation advertises: the worked example, the verification scripts, the
figure scripts, the pipeline runner, the data builders, following every `scripts/<name>.py`
reference transitively. A script is in the package because something advertised reaches it, not
because it happens to sit in the directory.

Some p82v1-named scripts are therefore included: the p84v2 runners call them, so they *are* the
implementations behind the reported results. Their names are kept rather than renamed, because
renaming a file that a frozen artifact's provenance string points to breaks the audit trail.

## What was deliberately excluded

| excluded | size | why | recoverable from |
|---|---|---|---|
| `data/structures/`, `data/raw/` | 917 MB | third-party structures; redistribution is unnecessary and the terms are not ours to grant | STCRDab and the PDB, `data/SOURCES.md`, `THIRD_PARTY_NOTICES.md`, `scripts/01_fetch_stcrdab.py` |
| `md/runs/` | 7.7 GB | MD trajectories; an open author decision, and too large to bundle | the included `.mdp` protocol, seeds and build manifest |
| `md/systems/*/charmm-gui*` | 3.3 GB | CHARMM-GUI build bundles | `md/systems/BUILD_MANIFEST.tsv` at charmm-gui.org |
| `results/checkpoints/` | 160 MB | trained weights from the **superseded 28-position candidate set**, under partitions that are not the frozen protocols. They are not the reported runs and reproduce nothing | nothing, and `results/review/FR4_CHECKPOINT_AUDIT.md` records why no usable checkpoint exists |
| `results/tables/_partial/p82v1/24b_run_all/` and `34b_peptide_aware/` | 72 MB | resume shards of a superseded generation. The **small** p82v1 partials are kept, because `55_stats_p84v2.py` validates its cluster-macro estimator against the recorded p82v1 values before applying it to p84v2 and fails without them | retraining, if ever wanted |
| `results/tables/_partial/p84v2/24b_run_all/` | 53 MB | resume shards of the pairwise runner; its *outputs* are included in full | retraining |
| `results/review/train_p84v2_logs/forced_backup_*/` | 94 MB | a backup taken during a re-run, superseded by the results themselves | n/a |
| `results/review/teim_determinism/{backup,runA}/` | 8 MB | superseded runs of the determinism investigation; the comparison tables and the report are included | n/a |
| `md/runs/`, `md/systems/*/charmm-gui*` | 11 GB | trajectories and build bundles, as above. The **cleaned input structures** `md/systems/*/prep/*.pdb` (2 MB) are kept, because Figure 7 reads them and they are the MD inputs | rebuild from `BUILD_MANIFEST.tsv` |
| `docs/manuscript/si_tables/superseded_20261102/` | small | SI tables from before the bootstrap correction; the corrected tables carry the superseded values in labelled columns | n/a |
| `__pycache__/`, `*.pyc`, `.ipynb_checkpoints/`, `queue.lock*` | small | machine-specific caches and runtime locks | n/a |
| the project's internal working notes and handoff documents at the repository root, `*_backup_*.md`, `LOG.md`, `PROJECT_STATUS.md`, `REVISED_PLAN.md`, the draft `*.md`/`*.docx` from earlier generations | n/a | internal working notes and session handoffs. They are the project's private process record, not part of the science, and several contain absolute local paths | the repository's own history |
| `ci5c00055.pdf`, `ci6c00337.pdf` | 14 MB | third-party published papers, kept locally for reference. Redistributing them is a copyright matter | their DOIs, in the reference list |
| the manuscript and SI `.docx`/`.pdf` | 17 MB | not published; distributing a manuscript is a publisher-rights question, not a reproducibility one | the **build assets** are included: both builders, `FIGURE_CAPTIONS.md` and the SI table generators rebuild both documents |
| `notebooks/`, `manuscript/`, `priors/`, `data/splits/` (v1) | n/a | earlier generations, superseded by `data/splits_v2/` and the current scripts | the repository's own history |

Nothing excluded is required by any command in `REPRODUCE.md`. **That is checked, not assumed.** The
staging script extracts every literal repository path that any shipped script reads and fails unless
each one is present in the package, so a file dropped by an exclusion rule is caught at staging
rather than discovered by a reader whose command will not run. Exactly four paths are allowed to be
absent, each with its reason recorded in the script: the third-party Struct-TCS identifier list
(referenced, not bundled), two output paths of a superseded figure generator, and a deliberately
nonexistent path in an MD smoke test that asserts a missing checkpoint is refused.

Four gaps were found this way and closed rather than documented away: `55_stats_p84v2.py` needs the
small p82v1 ablation partials for its estimator validation, Figure 7 needs the superseded
characterisation table whose SHA-256 it asserts against the run manifest and the four cleaned MD
input structures, and `67_si_control_tables.py` was crashing on a schema that had moved under it.

## What is in the package but is provenance rather than input

* `results/review/*.md`: the dated plans, the verification reports and the audit records, including
  the control plan committed before any downsampled model was trained, the comparator deviation
  table, the sweep-uncertainty reconciliation and the checkpoint audit. These are how a reader can
  tell what was decided before the results existed.
* `figures/previous_versions/` (Zenodo only): superseded figure versions, kept so the revision
  history in `figures/main/Figure_1.SOURCE.md` resolves to real files.
* `figures/source/fig1_approved_20261102/`: the approved Figure 1 exactly as the author supplied it,
  unmodified. It is the only shipped code that contains an absolute home path, deliberately: it is a
  historical record, and the **active** source `figures/source/fig1_arial_20261102/` is portable. The
  path check in `tests/test_frozen_artifacts.py` exempts that one directory by name and no other.

## Scientific provenance preserved

| question a reader may ask | where it is answered |
|---|---|
| What was fixed before the results existed, and when? | `results/review/CORRECTION_PLAN_P84V2.md`, `SIZE_MATCHED_CONTROL_PLAN.md`, `data/splits_v2/protocols.json` |
| How does each comparator differ from its publication? | `results/review/COMPARATOR_PROVENANCE.md` |
| Why do two interval estimates of the same quantity differ? | `results/review/SWEEP_UNCERTAINTY_RECONCILIATION.md` |
| Which figure came from which script and which table? | `figures/FIGURE_MANIFEST.md` |
| Which Figure 1 is in the manuscript, and what preceded it? | `figures/main/Figure_1.SOURCE.md` |
| What does each released array contain? | `DATA_DICTIONARY.md` |
| Why can one robustness check not be completed? | `results/review/FR4_CHECKPOINT_AUDIT.md` |
| What is third-party, and how is it obtained? | `THIRD_PARTY_NOTICES.md`, `data/SOURCES.md` |
