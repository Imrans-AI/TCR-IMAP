# Publication checklist

Everything below is **outstanding**. Nothing in this package has been pushed, uploaded, published or
deposited, no remote is configured, no release tag exists and no DOI has been minted or reserved.
This list is what turns the two local deliverables into a release.

## A. Author decisions — nothing can proceed without these

- [ ] **Licence.** `LICENSING_DECISION_REQUIRED.md`. Zenodo rejects an open deposit without one.
      Recommendation: MIT for the code, CC BY 4.0 for the derived data.
- [ ] **Repository destination.** Account or organisation, repository name, public or private at
      first push. No URL is written anywhere in this package until it exists.
- [ ] **ORCIDs** for both authors, or a decision to omit the field. None is recorded anywhere in
      this project, so none has been supplied.
- [ ] **Funding.** No funding statement exists in the project. Either supply grant identifiers for
      `ZENODO_METADATA.json` and a Funding section in the manuscript, or confirm there is none.
- [ ] **Copyright holder line**: the authors, Shanghai Jiao Tong University, or both.
- [ ] **Trajectory access.** The twelve 200 ns replicates (7.7 GB here, far more at full
      resolution) are not in the package. Decide: a second large Zenodo record, an institutional
      store, on request, or not available. **Until this is decided the manuscript must not promise
      any arrangement**, and at present it says only that the trajectories are available from the
      corresponding author on request — which is itself a commitment to honour.
- [ ] **Combined or split deposit.** This draft assumes one Zenodo record of type `dataset`. The
      alternative is two records, `software` and `dataset`, cross-linked.

## B. Before pushing the GitHub package

- [ ] Re-run `python release/packaging/stage_release.py` so both trees carry the final revision.
- [ ] `python tests/test_frozen_artifacts.py` → `8/8 checks passed`.
- [ ] `python scripts/58_verify_captions.py` → `ALL CAPTION NUMBERS VERIFIED`.
- [ ] `python scripts/worked_example.py` from a copy outside this repository.
- [ ] Add `LICENSE` (and `LICENSE-DATA` if the split is adopted).
- [ ] Fill `license` in `CITATION.cff`; add `repository-code` once the URL exists.
- [ ] Confirm nothing personal or machine-specific is in the tree: `AUDIT_PROVENANCE.md` lists what
      was excluded and why, and a check asserts that no shipped script contains an absolute home
      path.
- [ ] Create the release tag — the proposal is `v1.0.0` at acceptance, not now.

## C. Before depositing the Zenodo package

- [ ] Create the record, upload through the **bucket** API, then PUT the metadata, then publish.
      The legacy `/files` endpoint caps at 100 MB per file.
- [ ] Set `publication_date` to the actual deposit date.
- [ ] Set `license`.
- [ ] Resolve or delete every `REQUIRED-AUTHOR-DECISION` and `SET-ON-DEPOSIT` marker in
      `ZENODO_METADATA.json`. **A marker left in a published record is worse than an absent field.**
- [ ] Add the repository URL as an `isSupplementedBy` related identifier.
- [ ] Verify `sha256sum -c CHECKSUMS.sha256` passes on the uploaded archive after download.
- [ ] Record the minted DOI in `CITATION.cff`, `README.md` and the manuscript's availability section.

## D. Manuscript, once the records exist

- [ ] Replace the availability text in **both** places — Methods 2.12 and the Data and Code
      Availability section — with the repository URL and the archive DOI. Both currently state
      plainly that no deposit has been made and that no URL or identifier can be quoted; that is
      true today and must stay true until the records exist.
- [ ] Add the DOI to the Supporting Information statement if the journal wants it there.
- [ ] Settle the figure licence against the signed publisher agreement.
- [ ] Resolve the two items that are open for reasons unrelated to packaging: the **EPACT full
      author list** (reference 26) and the **GitHub handle** placeholder in the availability text.

## E. Standing constraints that a release must not break

- [ ] **Do not create or modify `md/analysis/RUN_AUTHORISED`.** The occupancy analysis is run-once
      by design within this project. The gate is deliberately not in the package; the analysis is
      fully documented so anyone else can repeat it, and no public-facing text says it must never be
      repeated.
- [ ] **Do not regenerate Figure 1.** It is an approved asset; the builder verifies its SHA-256 and
      aborts on a mismatch, and the superseded generator refuses to run without
      `--allow-superseded`.
- [ ] **Do not claim per-residue checkpoints or predictions exist.** They do not, and the limitation
      is recorded accurately in `results/review/FR4_CHECKPOINT_AUDIT.md`, the README and the release
      notes.
- [ ] **Do not regenerate the frozen protocols, the complex table or the contact labels.** They
      define the benchmark; changing one is version 2.0.0, not a patch.
