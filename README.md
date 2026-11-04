# TCR-IMAP: a dual-encoder cross-attention model for TCR–MHC platform contact map prediction

**Private review repository:** [Imrans-AI/TCR-IMAP](https://github.com/Imrans-AI/TCR-IMAP).
The code package was pushed privately on 3 November 2026. No public release or Zenodo deposit
exists, and no archive DOI has been assigned. See
`PUBLICATION_CHECKLIST.md` for what remains to be done, and `PACKAGE_INFO.json` for the exact source
revision this tree was staged from.

Accompanies the manuscript *TCR-IMAP: A Dual-Encoder Cross-Attention Model for TCR–MHC Platform
Contact Map Prediction*, by Muhammad Imran, Waqar Ali, Syed Mansoor Jan, Muhammad Rizwan, Adil
Farooq, Wanwan Huang and Dongqing Wei. The manuscript is being submitted for publication.

---

## What this is, and what it is not

**The contribution is the evaluation framework.** The claim is not that TCR-IMAP is a better
architecture; it is that how much of a sequence-based TCR–MHC contact predictor's apparent accuracy
belongs to the evaluation rather than to the model can be measured, and that this package contains
the frozen apparatus for measuring it: 174 curated structures, five partitioning protocols fixed
before the results were produced, contact labels and candidate masks keyed by IMGT number and
insertion code, within-fold frequency priors, paired per-complex metrics, a clustered bootstrap, and
every stored prediction the reported numbers were computed from.

Concretely, you get:

* **A characterized benchmark.** 174 TCR–peptide–MHC class I complexes, all with a human MHC class I
  heavy chain, collapsing to 69 TCR clusters at a similarity ratio of 0.70. The redundancy is
  measured, not assumed: 57.5 % of complexes have an exact TCR twin.
* **Five frozen partitioning protocols** – random, allele, temporal, tcr70 and a legacy
  intersection-selected partition retained from earlier work – with a declared primary metric,
  protocol and comparison. The legacy partition is **not** comparable to the others by strictness,
  because its held-out groups remain in training; it is reported for continuity, not as a rung on a
  scale. tcr70 is the primary protocol.
* **Baselines that are hard to beat for the right reason**: within-fold pair-frequency priors and
  nearest-neighbour label transfer, fitted on training folds only.
* **TCR-IMAP itself**, a dual-encoder cross-attention model that maps amino-acid identities to the
  complete 75 × 84 TCR-residue × MHC-residue contact map, plus a per-residue variant matching the
  output format of published residue-level predictors.
* **Comparator reimplementations** (Struct-TCS, Seq-TCS, TEIM-Res) with a provenance and deviation
  table recording every difference from the published specification.

### Supported tasks

| task | what it means here | reproducible from this package without a GPU |
|---|---|---|
| pairwise contact map | score every one of the 75 × 84 candidate residue pairs of a cohort complex | yes, from stored predictions |
| per-residue contact | per TCR position, does it contact any window residue | yes, from stored per-complex metrics |
| protocol comparison | the same methods under five frozen partitions, paired by complex | yes |
| identity-threshold sweep | accuracy as a function of the train/test similarity ceiling | yes |
| docking-angle regression | AngleNet under grouped evaluation | yes, from stored run records |
| MD occupancy | do high-scoring crystal-negative pairs form transient contacts in solution | yes, from the released occupancy tables |

### Limitations: read these before reusing anything

1. **There is no sequence-to-input builder.** `src/data/p84v2.py` exposes `T_for(pdbs)` and
   `M_for(pdbs)`, which look complexes up **by PDB identifier** in the frozen representation.
   Applying the model to a receptor outside this cohort needs an IMGT numbering step that this
   package documents but does not implement: ANARCI or equivalent for the TCR V domains, and
   per-domain alignment for the MHC G domain: never a fixed offset, which mismatches 66 of 126 α2
   positions. **Insertion codes must be preserved**: keying residues by IMGT number alone silently
   overwrites 12 residues at position 111, drops 63 at 112, and puts 1061A and 1072A in the base
   columns. That defect existed in an earlier build of this project and is the most likely way to
   misuse the representation.
2. **Scope of the evidence.** Paired αβ TCRs on human MHC class I, α1/α2 platform. Receptors are
   human in 165 of 174 complexes, mouse in 8, and one has a mixed chimpanzee/human α chain. MHC
   class II, γδ receptors and non-human MHC are out of scope.
3. **The predeclared primary comparison is unresolved at this sample size**: per-residue macro-AUPRC
   under tcr70, TCR-IMAP versus Struct-TCS, +0.0312 [−0.0049, +0.0650], p = 0.088. It is reported as
   it stands. No more favourable secondary comparison was promoted in its place.
4. **One robustness check could not be completed from stored artifacts.** The canonical-CDR-only
   (FR4) sensitivity covers 14 comparisons, the largest change being 0.0014, but **not** the
   predeclared per-residue primary: the per-residue benchmark writes no checkpoint at all and stores
   metrics rather than predictions, so that one comparison cannot be rescored without retraining.
   `results/review/FR4_CHECKPOINT_AUDIT.md` records the evidence. The gap is in one robustness
   check, not in the primary comparison, which was computed over its full predeclared candidate
   universe. **No historical per-residue checkpoint or prediction array exists, and none is claimed
   to.**
5. **The `double` protocol is legacy.** It selects the intersection of held-out cluster and allele
   groups but trains on all non-test complexes, so its held-out groups remain in training. It is
   kept, clearly labelled "intersection (legacy)", and is not part of any strictness ordering.
6. **Molecular-dynamics trajectories are not included**: see "Molecular dynamics" below.

---

## Quick start

```bash
conda env create -f env/environment.yml && conda activate tcr-imap   # or: pip install -r env/requirements.txt
python tests/test_frozen_artifacts.py       # integrity of the frozen artifacts, seconds
python scripts/worked_example.py            # the headline pairwise numbers, under a minute, CPU
```

`scripts/worked_example.py` checks the released arrays against the frozen labels, masks,
identifiers and cluster assignments, then recomputes the pairwise tcr70 comparisons and the identity
sweep, printing each computed value beside the value in the paper. CPU only, no network, no GPU.

If `python scripts/worked_example.py` reports missing arrays, this is the **GitHub tree** without the
frozen data. See `DATA.md`: the data lives in the Zenodo archive, and the two correspond by source
revision.

`REPRODUCE.md` lists every command that reproduces a table or figure from stored results, and gives
separate instructions and resource estimates for the optional retraining.

---

## Contents

| path | what it is |
|---|---|
| `src/` | library code: representation, metrics, bootstrap, priors, figure style |
| `scripts/` | the analysis and figure scripts, as a dependency closure of the advertised commands |
| `tests/` | integrity checks over the frozen artifacts; run them first |
| `configs/` | model and ablation configuration |
| `env/` | environment specification and the external tools that are not pip-installable |
| `data/splits_v2/protocols.json` | the five frozen partitioning protocols, never regenerated |
| `data/processed/complex_metadata.tsv` | the 174-complex table with clu70/80/90, never regenerated |
| `data/processed/pair_maps_p84v2.npz` | contact labels and masks, 174 × 75 × 84, with the row and column manifests |
| `data/interim/stcrdab_index.tsv` | the STCRDab index rows for the cohort |
| `results/tables/` | stored predictions, per-complex metrics, aggregates and the MD output tables |
| `results/review/` | the control, sweep, sensitivity and comparator analyses with their inputs and dated plans |
| `figures/` | the final main and supporting figures, the manifest, and the approved Figure 1 source |
| `docs/` | figure captions, SI tables and their builders, the manuscript and SI build assets |
| `md/` | the MD selection criteria, protocol, run scripts, analysis code and its own tests |
| `DATA_DICTIONARY.md` | every released array and table, field by field |
| `FILE_MANIFEST.tsv`, `CHECKSUMS.sha256` | size and SHA-256 of every file in this package |
| `AUDIT_PROVENANCE.md` | what was included, what was left out and why |
| `PACKAGE_INFO.json` | version, staging date, source revision, Figure 1 hash, deposit status |

---

## Molecular dynamics

`md/MD_SYSTEM_SELECTION_CRITERIA.md` holds the predeclared criteria, the system selection, the
protocol and the analysis definitions. `md/protocol/*.mdp` are the GROMACS parameter files,
`md/systems/BUILD_MANIFEST.tsv` records how the four systems were built, `md/analysis/` holds the
analysis code and its tests, and `results/tables/md_*.csv` hold the twelve output tables the single
authorised run produced. An independent party can therefore repeat the occupancy analysis.

Within this project the occupancy analysis is run-once by design, gated by a local authorisation
file that is deliberately not included. **That gate is a project-internal safeguard against
re-running a predeclared analysis after seeing its result, not a restriction on anyone else, and
nothing here says the analysis must never be repeated.**

**Trajectories and CHARMM-GUI bundles are not in this package.** The twelve 200 ns replicates are
7.7 GB as stored here and far larger with full-resolution coordinates, and the four CHARMM-GUI
bundles are 3.3 GB. Everything needed to regenerate them is included: the structures are public, the
build manifest records the CHARMM-GUI parameters, and the `.mdp` files carry the protocol and seeds.
The MD trajectories are available from the corresponding author on reasonable request.

---

## Third-party material

See `THIRD_PARTY_NOTICES.md`. The full source-structure collection is retrieved from STCRDab
and the Protein Data Bank. Four cleaned coordinate files (3vxs, 5yxn, 7n1f and 8eo8) are included
as MD inputs; the remaining coordinate files are referenced through retrieval instructions. The comparator
reimplementations are our own code written from the published descriptions; no third-party source
repository is bundled, and `results/review/COMPARATOR_PROVENANCE.md` records every deviation from
each published specification.

## Licence

Original code is licensed under **MIT** (`LICENSE`). Original derived data and documentation
are licensed under **CC BY 4.0** (`LICENSE-DATA`). Copyright 2026 Muhammad Imran and Dongqing Wei.
See `LICENSING_DECISION_REQUIRED.md` for scope and exclusions: third-party material retains its
own terms. The final figures are all original, generated by scripts in this repository, and are
released under CC BY 4.0 with the other original derived material.

## Citing

`CITATION.cff` carries the authors, title, repository URL and applicable licences. No DOI or
public release date is supplied because neither exists yet.

## Funding and acknowledgments

D.W. is supported by the National Key R&D Program of China (2024YFA1306400 and 2023YFE0199200), the National Natural Science Foundation of China (32570803 and 32030063), and the Intergovernmental International Scientific and Technological Innovation Cooperation Program. The computations were performed in part at Pengcheng Laboratory and the Center for High-Performance Computing, Shanghai Jiao Tong University.
