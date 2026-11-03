# Licensing — an open author decision

**No licence is assigned, and no licence file is shipped.** Nothing in this repository records a
licensing decision: there is no `LICENSE`, no `COPYING`, no licence header in any source file and no
statement in `CLAUDE.md`, `GUIDELINES.md` or the manuscript. Assigning one is the author's call, and
guessing would be worse than leaving it open, so this file states the recommendation and stops.

**Until a licence is chosen, default copyright applies**: all rights reserved, no permission to
reuse, and the package is for author review only. Zenodo additionally **requires** a `license` value
whenever `access_right` is `open`, so a deposit cannot be published without this decision.

## Recommendation

A split licence, which is the normal arrangement for a paper with code and data, and the one the
journals and the archive expect:

| part | recommended | why |
|---|---|---|
| **code** — `src/`, `scripts/`, `tests/`, `md/analysis/`, `figures/source/`, the manuscript builders | **MIT** | permissive, one page, imposes nothing on reusers, and matches the comparator reimplementations' purpose of being inspected and adapted. BSD-3-Clause is equivalent in effect if an institution prefers it |
| **derived data** — contact labels, masks, the cohort table, the frozen protocols, stored predictions, per-complex metrics, MD output tables | **CC BY 4.0** | the frozen protocols and the characterized benchmark are the contribution others are most likely to reuse, and attribution is the only condition worth imposing. CC0 is the alternative if the aim is maximum reuse with no conditions; Zenodo defaults datasets to CC0 |
| **figures** — `figures/main/`, `figures/supplementary/` | follow whatever the journal's copyright-transfer agreement allows, typically **CC BY 4.0** | this is the one part a publisher agreement can constrain, so it must be settled after the agreement is signed, not before |

Avoid a non-commercial (`-NC`) or no-derivatives (`-ND`) clause on the data. Both make the benchmark
unusable for the reuse the paper argues for, and `-NC` is widely read as blocking industry research.

## What the decision affects

1. `LICENSE` and, if the split is adopted, `LICENSE-DATA` at the package root.
2. The `license` field in `ZENODO_METADATA.json`, which is currently
   `REQUIRED-AUTHOR-DECISION`. Zenodo expects an identifier from its `/api/licenses` list, for
   example `mit` or `cc-by-4.0`.
3. The `license` field in `CITATION.cff`, currently absent for the same reason.
4. A one-line "Licence" section in `README.md`, which currently says the decision is open.
5. Whether a licence header is added to source files. Not required by MIT; a `SPDX-License-Identifier`
   line is a cheap convention if wanted.

## What the decision does not affect

The third-party material in `THIRD_PARTY_NOTICES.md`. The structures stay under their own terms and
are referenced rather than redistributed, so no licence chosen here can or does apply to them.

## Required of the author

- [ ] Confirm the code licence (recommended **MIT**).
- [ ] Confirm the data licence (recommended **CC BY 4.0**).
- [ ] Confirm that the institution and the funder impose no conflicting requirement.
- [ ] Confirm the figure licence **after** the publisher agreement is in place.
- [ ] State the copyright holder line: the authors, or Shanghai Jiao Tong University, or both.
