# Module 4 → JCIM: standing guidelines

*From handwritten notes, 19 Aug 2026. These are directives for how the paper must be
strengthened before submission — read alongside the project's internal notes.*

---

## 1. Make the architecture strong

The model must not read as an incidental analysis tool. Concretely:

- Run **on the competitor's data or splits** — i.e. reproduce the SG-TPMI / Struct-TCS
  STRUCT-CS setting (their 94 complexes, their partition) and show GRIP performs at least
  comparably there, *in addition* to our stricter protocols. This pre-empts "you only win
  on your own benchmark."
- Current honest position: GRIP > Struct-TCS on tcr70 (+0.0273, p=0.001) but parity elsewhere;
  ablations show much of the signal is positional (`shuffle_tcr` n.s., `no_cross` n.s.).
  Strengthening the architecture means finding a component whose removal *significantly* hurts.

**Open task:** rebuild the STRUCT-CS-equivalent evaluation from their PDB list
(`data/raw/sgtpmi/data/STRUCT-CS/PDBids.txt`) and report side by side.

## 2. Confirm the splitting strategy

Splits are frozen in `data/splits_v2/protocols.json` — but the *justification* must be airtight:

- tcr70 is primary (0% of test complexes with a >90%-identical training neighbour, 0 duplicate-family leakage)
- Report the similarity audit table as a figure, not a footnote
- Sensitivity across 70 / 80 / 90% clustering thresholds must be shown, not just the best one
- Double-holdout has *higher* measured leakage (0.422) than tcr70 (0.000) — say so plainly rather
  than implying "double" is strictest
- Every prior/frequency baseline must be fit **within training folds only** — re-verify before submission

## 3. Match JCIM-level figure standards

Reference: SG-TPMI (JCIM 2025) and HLABrew (JCIM 2026) — dense, multi-panel, two-column figures.

**Structures must be emphasised.** Currently missing entirely:
- PyMOL renderings of representative TCR–pMHC complexes
- The recovered α1/α2 helix footprint projected onto an actual structure surface
- Contact maps shown *on* the structure, not only as heatmaps
- Representative success and failure cases with predicted vs observed contacts side by side

**Contacts must be emphasised.** Contact geometry is the paper's subject; the figures should
make it visible rather than abstract.

Planned figure set:
1. Dataset construction and evaluation design (schematic; BioRender)
2. Dataset redundancy and train–test similarity per protocol
3. Similarity-threshold curve *(built — needs polish)*
4. Benchmark vs baselines across protocols with bootstrap CIs *(built — needs polish)*
5. Contact stereotypy: heatmap + **structure projection** + entropy with matched resampling *(partial)*
6. External/temporal test and failure analysis

Figures should be assembled in Illustrator/Inkscape after plotting; matplotlib defaults will not
reach this standard unaided.

---

## Immediate next actions

- [ ] Reproduce the STRUCT-CS split and benchmark GRIP on the competitor's own data
- [ ] Stereotypy entropy with matched resampling + bootstrap CIs (values in Fig 5b are still hardcoded)
- [ ] PyMOL structure panels for Figures 1, 5 and 6
- [ ] Peptide-aware variant + peptide-shuffled control
- [ ] Archive the SG-TPMI commit and publish the reimplementation patch
