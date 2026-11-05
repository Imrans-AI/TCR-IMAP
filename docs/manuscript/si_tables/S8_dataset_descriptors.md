# Table S8. Dataset descriptors

**Table S8.** Descriptive statistics of the 174-complex benchmark that are quoted in the text but were previously
recorded only in review documents, each with the computation that produces it.

| Quantity | Value | Computation | Source |
|---|---|---|---|
| Exact TCR twins | 100 of 174 complexes (57.5 %) | Complexes whose nearest neighbor (another complex) has a Levenshtein similarity ratio ≥ 0.999 on the clustering strings. | `tcr_identity_matrix.npy`; `identity_bound/handoff6_statements.tsv` (script `54b_dupfamily_and_handoff6.py`) |
| Twin pairs differing in peptide | 123 of 163 twin pairs (75.5 %) | Among all complex pairs with similarity ≥ 0.999, those whose `peptide` field differs. | `tcr_identity_matrix.npy`, `complex_metadata.tsv`; `handoff6_statements.tsv` (54b) |
| Deduplication is insufficient | 42 complexes removed → 132 remain; median nearest-neighbor similarity 0.967 | Remove every complex in a duplicate family with more than one member (the stored `dup_family`: identical CDR3 pair, allele and peptide); for each remaining complex, its maximum similarity to any other remaining complex; median over complexes. | `complex_metadata.tsv` (`dup_family`), `tcr_identity_matrix.npy`; `handoff6_statements.tsv` (54b) |
| Helix share of G-domain contacts | 4,674 of 5,294 contacts (88.3 %) | TCR–MHC contacts (minimum heavy-atom distance ≤ 4.5 Å) over the 182-position α1/α2 G-domain, summed over the six loops, on positions whose DSSP majority class is α1 or α2 helix. | `gdomain182_region_loop_summary_p84v2.csv` (script `53_gdomain182_descriptive_p84v2.py`); `REBUILD_P84V2_REPORT.md` step 4 |
| β-sheet floor | 0 contacts in 595,952 resolved pairs | Same table, positions classified β-sheet floor (56 positions). | `gdomain182_region_loop_summary_p84v2.csv` |

**Notes.**
- Similarity: Levenshtein similarity ratio on the clustering strings (CDR1/CDR2/CDR3 slot segments of both chains, IMGT
  27–38, 56–65 and 105–124, FR4 118–124 included, one residue per IMGT position), the strings that define clu70 and the
  frozen splits. Under the other string definitions of `IDENTITY_BOUND_REPORT.md` the twin and peptide statistics are
  unchanged and the deduplication median is 0.968 (full CDR-plus-framework sequences) or 0.958 (canonical CDRs alone).
- The deduplication uses the stored, allele-dependent `dup_family` (CDR3 pair, allele, peptide), as the original
  statement did; its grouping over 108 missing alleles is a known limitation of the frozen metadata (Table S4).
- The helix share is computed over the 182-position G-domain (IMGT 1–90 and 1001–1090 plus 1061A and 1072A), of whose
  5,294 contacts the 84-column window holds 5,281. Region classes are DSSP majority calls (helix H/G/I, sheet E/B,
  ≥ 30 resolving structures, majority > 0.5).
- Every value was recomputed by `build_si_tables.py` from the committed inputs and asserted equal to the recorded value.

**Sources.** Listed per row. Recorded values: `results/review/identity_bound/handoff6_statements.tsv` (definition 1),
`results/review/REBUILD_P84V2_REPORT.md` step 4; the remaining recorded values are those asserted during the build.
