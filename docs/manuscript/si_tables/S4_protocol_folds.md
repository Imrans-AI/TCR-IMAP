# Table S4. Partitioning protocols: fold composition

**Table S4.** Composition of every fold of the five frozen partitioning protocols: complexes and clu70 clusters on each
side, the train/test TCR similarity, and duplicate families split across the fold boundary.

| Protocol | Fold | Evaluated | n test | n training | Test clusters | Training clusters | Clusters on both sides | Max train/test similarity | Mean of per-test maximum | Test ≥ 0.90 (%) | Dup. families split (CDR3 pair + peptide) | Dup. families split (+ allele, stored) | Per-residue models: training / validation (validation clusters; % held out) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| random | 0 | yes | 35 | 139 | 25 | 64 | 20 | 1.000 | 0.911 | 71.4 | 9 | 8 | 116 / 23 (11; 16.5 %) |
|  | 1 | yes | 35 | 139 | 25 | 58 | 14 | 1.000 | 0.829 | 48.6 | 6 | 5 | 118 / 21 (13; 15.1 %) |
|  | 2 | yes | 35 | 139 | 26 | 61 | 18 | 1.000 | 0.880 | 60.0 | 7 | 5 | 118 / 21 (5; 15.1 %) |
|  | 3 | yes | 35 | 139 | 26 | 57 | 14 | 1.000 | 0.860 | 60.0 | 9 | 9 | 111 / 28 (10; 20.1 %) |
|  | 4 | yes | 34 | 140 | 20 | 64 | 15 | 1.000 | 0.901 | 61.8 | 6 | 5 | 117 / 23 (9; 16.4 %) |
| intersection (legacy) | 0 | yes | 16 | 158 | 10 | 61 | 2 | 0.850 | 0.671 | 0.0 | 0 | 0 | 134 / 24 (7; 15.2 %) |
|  | 1 | yes | 5 | 169 | 4 | 67 | 2 | 1.000 | 0.851 | 60.0 | 0 | 0 | 142 / 27 (17; 16.0 %) |
|  | 2 | yes | 9 | 165 | 8 | 67 | 6 | 1.000 | 0.902 | 66.7 | 2 | 0 | 138 / 27 (17; 16.4 %) |
| allele | 0 | yes | 26 | 148 | 16 | 62 | 9 | 1.000 | 0.837 | 50.0 | 3 | 0 | 122 / 26 (11; 17.6 %) |
|  | 1 | yes | 20 | 154 | 11 | 62 | 4 | 1.000 | 0.768 | 25.0 | 1 | 0 | 123 / 31 (13; 20.1 %) |
|  | 2 | yes | 15 | 159 | 10 | 62 | 3 | 1.000 | 0.714 | 13.3 | 0 | 0 | 130 / 29 (8; 18.2 %) |
|  | 3 | **no (n < 3)** | 2 | 172 | 2 | 68 | 1 | 1.000 | 0.828 | 50.0 | 1 | 0 | — |
|  | 4 | yes | 3 | 171 | 3 | 67 | 1 | 1.000 | 0.783 | 33.3 | 1 | 0 | 145 / 26 (8; 15.2 %) |
| temporal | 0 | yes | 60 | 114 | 32 | 45 | 8 | 1.000 | 0.688 | 5.0 | 0 | 0 | 91 / 23 (10; 20.2 %) |
| tcr70 | 0 | yes | 39 | 135 | 14 | 55 | 0 | 0.689 | 0.634 | 0.0 | 0 | 0 | 113 / 22 (7; 16.3 %) |
|  | 1 | yes | 33 | 141 | 14 | 55 | 0 | 0.696 | 0.664 | 0.0 | 0 | 0 | 119 / 22 (10; 15.6 %) |
|  | 2 | yes | 25 | 149 | 14 | 55 | 0 | 0.696 | 0.667 | 0.0 | 0 | 0 | 126 / 23 (5; 15.4 %) |
|  | 3 | yes | 32 | 142 | 14 | 55 | 0 | 0.695 | 0.663 | 0.0 | 0 | 0 | 120 / 22 (10; 15.5 %) |
|  | 4 | yes | 45 | 129 | 13 | 56 | 0 | 0.696 | 0.672 | 0.0 | 0 | 0 | 104 / 25 (13; 19.4 %) |

Complexes evaluated: random 174, intersection (legacy) 30, allele 64, temporal 60,
tcr70 174 of 174.

**Notes.**
- n training is every complex outside the test fold (174 − n test). The pairwise models and the frequency priors train on
  all of them. The per-residue models (TCR-IMAP per-residue, Struct-TCS, Seq-TCS) hold out whole clu70 clusters for
  validation until at least 15 % of those complexes are held out; the last column gives that split for evaluated folds.
- Similarity is the Levenshtein similarity ratio (`Levenshtein.ratio`), an edit-distance similarity, not gap-free alignment
  identity, computed on the clustering strings: the six CDR1/CDR2/CDR3 slot segments of both chains, IMGT 27–38, 56–65 and
  105–124 (the last including the J-region framework positions 118–124), one residue per IMGT position. "Max train/test
  similarity" is the largest, over the fold's test complexes, of each test complex's maximum similarity to any training
  complex. Under other string definitions the tcr70 maximum is 0.707 (full CDR-plus-framework sequences; 12 test
  complexes between 0.70 and 0.71, none above 0.80) and 0.682 (canonical CDR sequences alone); see Table S8 and
  `IDENTITY_BOUND_REPORT.md`.
- **Allele and intersection coverage limitation.** The allele and intersection protocols cover only 64 and 30 of the 174 complexes.
  Allele annotation is regex-parsed from the free-text `compound` field of the STCRDab summary, so 108 complexes have
  none. In `scripts/23_make_protocols.py` (lines 26 and 33) the allele groups for the folds were formed with
  `fillna("NA")`, but fold membership was then tested with `.isin()` on the original column, where those 108 complexes
  are still missing; they therefore enter no allele or intersection test fold. Allele fold 3 holds only two complexes and is
  dropped by the three-complex minimum used in every training script. The frozen protocol file was not regenerated.
- Duplicate families: identical CDR3α–CDR3β pair and peptide (allele-independent; the definition reported) and, as
  stored in `complex_metadata.tsv`, identical CDR3 pair, allele and peptide (allele-dependent; grouped over 108 missing
  alleles, so it both over- and under-merges). tcr70 is 0 in every fold under both.
- Fold-level duplicate counts include allele fold 3; the mean over the four evaluated allele folds is 1.25 (Figure 2E),
  and 1.2 over all five.

**Sources.** `data/splits_v2/protocols.json` (SHA-256 verified), `data/processed/complex_metadata.tsv` (clu70; SHA-256
verified), `results/tables/tcr_identity_matrix.npy` (similarity), checked against `results/tables/leakage_audit_dual.tsv`
(n test, mean of per-test maximum; duplicate families taken from its `leak_conservative` and `leak_stored` columns) and
`results/tables/percomplex_maxsim.tsv` (per-fold maximum); per-residue validation split from
`results/review/validation_fix/B2_split_dryrun.log`. Computed by `build_si_tables.py` (arithmetic only).
