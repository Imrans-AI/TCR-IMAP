# Building the v3 documents

```bash
SP=/tmp/claude-1000/<session>/scratchpad          # the directory holding node_modules/docx
cd docs/manuscript
NODE_PATH=$SP/node_modules node build_manuscript_docx_v3.js
NODE_PATH=$SP/node_modules node build_si_docx_v3.js
```

## The docx version matters — this is a trap

**Use `docx` 9.15.0.** Version **9.7.1** silently breaks every equation.

`docs/manuscript/` has no `node_modules` on its own path, so Node falls back to `NODE_PATH`. Two
copies exist on this machine:

| location | version | result |
|---|---|---|
| the session scratchpad (`$SP/node_modules`) | **9.15.0** | correct — this built v2 and v3 |
| `~/PhD_data/TCR_Module_1/node_modules` | 9.7.1 | **breaks all 13 equations** |

Under 9.7.1 the OMML loses `xml:space="preserve"` on every math run and omits the empty
`<m:sup/>` that `<m:nary>` requires. The file still builds and the *text* of every equation is
correct, so no error is raised — but LibreOffice then renders the summations as `□` placeholders.
Caught by converting to PDF and reading the page; it is invisible in the XML text extraction.

**Check after any build:**

```python
import zipfile, re
x = zipfile.ZipFile("TCR-IMAP_manuscript_v3.docx").read("word/document.xml").decode("utf8")
e = re.findall(r"<m:oMath>.*?</m:oMath>", x, re.S)
assert len(e) == 13
assert sum('xml:space="preserve"' in q for q in e) == 13   # 0 under 9.7.1
assert sum('<m:sup/>' in q for q in e) == 5                # 0 under 9.7.1
```

Since the scratchpad is session-scoped, `npm install docx@9.15.0` into a stable directory before
the next build and point `NODE_PATH` at it.

## Title

Both builders carry the author-selected title **“TCR-IMAP: Homology-Controlled Evaluation of
TCR–MHC Contact Prediction”** (selected 2 November 2026), in the title block and in the DOCX
metadata. v1 and v2 keep their own titles and are not rebuilt. The sequence-only nature of the
model is carried by the abstract (“from sequence alone”, “sequence-based”) and by Methods 2.3 and
2.5 (“amino-acid identities only; no coordinates”); the shorter title does not widen the
demonstrated scope, which stays human MHC class I, α1/α2 platform, paired αβ TCRs.

## What v3 contains

Forked from the v2 builders at commit `202b61b`, so **v2 remains reproducible byte for byte** —
`build_manuscript_docx_v2.js` and `TCR-IMAP_manuscript_v2.docx` are untouched. Verified: rebuilding
the v2 builder reproduces the committed v2 equations exactly, and v3 differs from v2 in exactly one
equation (Equation 2, `d_L` → `d_I`).

Changes are listed in `results/review/V3_PROPOSED_TEXT.md`, with the evidence for each in
`results/review/V3_PART_A_VERIFICATION.md`, `V3_PART_C_VALIDATION_SPLIT.md`,
`COMPARATOR_PROVENANCE.md` and `SIZE_MATCHED_CONTROL_RESULT.md`.

## Figure 1 is an approved asset, and the build refuses to ship anything else

The builder embeds **`figures/main/Figure_1.png`**, the author-approved figure installed on
2 November 2026 and revised the same day to Arial typography (6,480 × 4,230 px, 7.2 × 4.7 in at
900 dpi, sha256 `0684c8647a0b158c86af9229f2344f1eec4f8ae71b08e83c12cda0a26364d9a1`; the superseded
DejaVu Sans version is `a447b2a8…b45bbd0`, preserved under `figures/previous_versions/`). It is scaled uniformly to
the 6.5 in text width (6.5000 × 4.2431 in, ≈ 997 dpi effective) with `exactAspect`, so the
proportions are preserved exactly and no pixels are discarded. Provenance, the editable source
and the restore copies are in `figures/main/Figure_1.SOURCE.md`,
`figures/source/fig1_arial_20261102/` (active) and `figures/source/fig1_approved_20261102/`
(as supplied).

`build_manuscript_docx_v3.js` hashes the file on every build and **throws** if it is not the
approved figure, so a regenerated or substituted Figure 1 cannot reach the DOCX silently.
`scripts/70_fig1_workflow.py`, which produced the superseded script-drawn Figure 1, is kept for
the audit trail but refuses to run without `--allow-superseded` and writes only to
`figures/previous_versions/`.

## Release decisions, 3 November 2026

The private repository is https://github.com/Imrans-AI/TCR-IMAP. The EPACT reference already carries all nine
authors; the former checklist entry was stale. The author approved MIT for original code and
CC BY 4.0 for original derived data/documentation, authors as copyright holders, the funding
statement in `release/packaging/docs/FUNDING.md`, and trajectories on reasonable request.
These changes are in the v3 builder and DOCX. The existing DOCX was edited directly at the OOXML
paragraph level to preserve all equation and media parts, then rendered for verification.
Figure assets have not been relicensed. No public release, tag or Zenodo deposit exists.

Validation of the 3 November edit: 37 pages; funding on page 34; only pages 15, 16 and 34
differ pixel-for-pixel from the prior DOCX rendered with the same runtime. All 13 equations,
all 63 citation targets and all embedded media remain unchanged. The explicit blank paragraph
carrying the References page break was replaced with pageBreakBefore on its heading to avoid
an empty overflow page after adding funding. SI and all historical v2 files were untouched.
