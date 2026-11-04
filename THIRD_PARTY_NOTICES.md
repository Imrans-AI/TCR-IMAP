# Third-party material

This package contains original code and derived data, plus the four cleaned coordinate files
identified below. Third-party material retains its own terms.

## Structures: source collection and included MD inputs

| material | source | how to retrieve |
|---|---|---|
| 174 IMGT-numbered TCR–peptide–MHC class I structures | **STCRDab**, Leem *et al.*, *Nucleic Acids Res.* 2018, 46, D406–D412, retrieved 15 August 2026 | https://opig.stats.ox.ac.uk/webapps/stcrdab/ ; `data/interim/stcrdab_index.tsv` carries the exact summary rows, and `scripts/01_fetch_stcrdab.py` performs the download |
| the underlying crystallographic entries | **Protein Data Bank** | https://www.rcsb.org ; the identifiers are the `pdb` column of `data/processed/complex_metadata.tsv` |

`data/SOURCES.md` records the retrieval date and the selection criteria. The **derived** data in
this package: contact labels, masks, residue graphs, the aligned representation, the cohort table
and the clusters: is our own work computed from those public structures.

Four cleaned coordinate files are included for the MD inputs and Figure 7:

- `md/systems/3vxs/prep/3vxs_clean.pdb`: PDB 3VXS
- `md/systems/5yxn/prep/5yxn_clean.pdb`: PDB 5YXN
- `md/systems/7n1f/prep/7n1f_clean.pdb`: PDB 7N1F
- `md/systems/8eo8/prep/8eo8_clean.pdb`: PDB 8EO8

These are processed copies of source coordinates, not newly determined structures. Preparation
is recorded in `md/systems/BUILD_MANIFEST.tsv`. The package's MIT and CC BY 4.0 grants do not
relicense third-party coordinates. The rest of the structure collection is not bundled.

One structure, `1ao7`, is rendered in Figure 1 panel A. The rendering is ours
(`figures/source/fig1_arial_20261102/render_structure.py`); the coordinates are the public PDB entry
and are not shipped.

## Comparator methods: reimplemented, not redistributed

The Struct-TCS, Seq-TCS and TEIM-Res comparators in this package are **our own code, written from
the published descriptions**. No original author's repository, weights or data file is included.

| comparator | published source |
|---|---|
| Struct-TCS, Seq-TCS | Gao *et al.*, *J. Chem. Inf. Model.* 2025, 65, 4674–4686, DOI 10.1021/acs.jcim.5c00055 |
| TEIM-Res | as cited in the manuscript reference list |

**One third-party file is referenced rather than bundled.** `scripts/31b_structcs_setting_p82v1.py`
and its pre-p82v1 predecessor read
`data/raw/sgtpmi/data/STRUCT-CS/PDBids.txt`: the list of PDB identifiers in the comparator's own
test set, taken from the released repository `gaol00034/pan-pMHC` (local checkout at commit
`c118be5`). It is 475 bytes and is **not** in this package, because permission to redistribute it is
unclear. Retrieve it from that repository if you want to re-run the comparator's own single-split
design; the **results** of that design, `results/tables/structcs_setting_p84v2_{4p5A,6A}.csv`, are
included, so nothing in the paper depends on having it. Which complexes the subset contains is also
described in `results/review/COMPARATOR_PROVENANCE.md`.

`results/review/COMPARATOR_PROVENANCE.md` is the deviation table: published versus implemented task
and inputs, residue features, graph construction, architecture, training and validation, the
execution fixes that were needed, the deliberate adaptations, and the components that were not
available at all. **This is an architecture comparison under common inputs, not a reproduction of
published performance**, and the manuscript says so. Where our implementation departs from a paper's
text or from its released code, the departure is ours and is recorded there, not attributed to the
original authors.

## Software

The dependencies in `env/environment.yml` are installed from PyPI and conda-forge under their own
licences and are not redistributed here. The external tools in `env/EXTERNAL_TOOLS.md`, GROMACS,
CHARMM-GUI, DSSP, PyMOL: are obtained from their own distributors under their own terms.

## Fonts

Figures embed subsets of **Arial** (`ArialMT`, `Arial-BoldMT`) through matplotlib's standard PDF
font subsetting, as published figures normally do. No font file is redistributed as a font.

## Licence of this package

MIT applies to original code and CC BY 4.0 to original derived data and documentation.
See `LICENSING_DECISION_REQUIRED.md` for the scope and exclusions. These grants do not extend
to third-party material or to the final figure assets.
