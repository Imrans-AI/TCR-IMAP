# Building the v4 documents

```bash
R=$PWD                                    # repository root
cd docs/manuscript
NODE_PATH=$R/.build_tools/node_modules node build_manuscript_docx_v4.js
NODE_PATH=$R/.build_tools/node_modules node build_si_docx_v4.js
```

`.build_tools/` holds the pinned `docx` package and is not part of this repository; create it with
`npm install --prefix .build_tools docx@9.15.0`.

The builders are the authority. **Never edit the `.docx` by hand** — any such edit is lost at the
next build, and see the LibreOffice warning below.

## The docx version matters — this is a trap

**Use `docx` 9.15.0.** Version **9.7.1** silently breaks every equation.

Under 9.7.1 the OMML loses `xml:space="preserve"` on every math run and omits the empty `<m:sup/>`
that `<m:nary>` requires. The file still builds and the *text* of every equation is correct, so no
error is raised — but the summations then render as `□` placeholders. It is invisible in an XML
text extraction; it shows up on the page.

## Never save the `.docx` from LibreOffice

Open it read-only. A LibreOffice save does not destroy the equations — all 13 survive as OMML — but
it rewrites them: intentional spacing inside the math is stripped, the hyphen in `macro-AUPRC`
becomes a minus sign (U+2212), and `1/B` is re-parsed into a fraction. It also renames the media
parts to `image1.png`…, losing the content-hashed names, and drops `docProps/custom.xml`,
`word/footnotes.xml` and `word/endnotes.xml`. Recover with `git checkout --`, or rebuild.

Converting to PDF with `soffice --headless --convert-to pdf` is safe: it does not modify the input.
Do it from a copy anyway, and check the `.docx` checksum is unchanged.

## Verify after every build

```python
import zipfile, re
z = zipfile.ZipFile("TCR-IMAP_manuscript_v_4.docx")
x = z.read("word/document.xml").decode("utf8")
assert len(re.findall(r"<m:oMath[ >]", x)) == 13      # allow attributes: see below
assert not any(n.startswith("word/media/image") for n in z.namelist())
```

**Count the tag allowing attributes.** LibreOffice writes `<m:oMath xmlns:m="…">`, so a grep for the
literal `<m:oMath>` returns 0 on a LibreOffice-saved file whose equations are intact. That false
negative has caused three separate false alarms in this project.

Expected for v4: **13 equations, 63 references numbered 1…63 contiguously, 41 pages** for the
manuscript and **31 pages** for the Supporting Information.

## Generated inputs

`build_manuscript_docx_v4.js` reads `docs/manuscript/results_tables_v4.json` for Tables 1 and 2 and
verifies a SHA-256 over each source file that produced it. If a source changes, the build stops and
tells you to re-run:

```bash
python scripts/71_results_tables_v4.py
```

`scripts/72_references_acs.py` records the ACS reference formatting and the PubMed/publisher
verification behind it.
