#!/usr/bin/env python3
"""Integrity checks over the frozen artifacts that every reported number depends on.

Runs with pytest or on its own:

    pytest -q tests/
    python tests/test_frozen_artifacts.py

Checks that need the frozen arrays are skipped, with a clear message, when only the code tree is
present (a GitHub checkout without the Zenodo data). Nothing here trains, scores or writes.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/pair_maps_p84v2.npz"
PAIRMAPS = ROOT / "results/tables/pairmaps_tcr70_p84v2.npz"


class _Skip(Exception):
    """Raised when a check needs an artifact that is only in the Zenodo archive."""


def _skip(msg):
    """Skip under pytest, or signal the standalone runner to record a skip."""
    try:
        import pytest
    except ImportError:
        raise _Skip(msg)
    pytest.skip(msg)


def test_metadata_cohort():
    """174 complexes, 69 clu70 clusters, unique PDB identifiers."""
    md = pd.read_csv(ROOT / "data/processed/complex_metadata.tsv", sep="\t")
    assert len(md) == 174, len(md)
    assert md.pdb.is_unique
    assert md.clu70.nunique() == 69
    assert md.clu80.nunique() == 82 and md.clu90.nunique() == 89


def test_protocols_frozen():
    """Five protocols; tcr70 evaluates all 174 complexes and no fold is below the minimum."""
    pr = json.loads((ROOT / "data/splits_v2/protocols.json").read_text())
    assert set(pr) == {"random", "double", "allele", "temporal", "tcr70"}, sorted(pr)
    evaluated = {k: sum(len(f) for f in v if len(f) >= 3) for k, v in pr.items()}
    assert evaluated["tcr70"] == 174, evaluated
    assert evaluated == {"random": 174, "double": 30, "allele": 64,
                         "temporal": 60, "tcr70": 174}, evaluated
    md = pd.read_csv(ROOT / "data/processed/complex_metadata.tsv", sep="\t")
    for name, folds in pr.items():
        flat = [p for f in folds for p in f]
        assert len(flat) == len(set(flat)), f"{name}: a complex appears in two folds"
        assert set(flat) <= set(md.pdb), f"{name}: unknown identifiers"


def test_tcr70_has_no_cluster_overlap():
    """The primary protocol holds out whole clu70 clusters."""
    pr = json.loads((ROOT / "data/splits_v2/protocols.json").read_text())
    md = pd.read_csv(ROOT / "data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
    folds = pr["tcr70"]
    for i, test in enumerate(folds):
        train = [p for j, f in enumerate(folds) if j != i for p in f]
        assert not (set(md.loc[test].clu70) & set(md.loc[train].clu70)), f"fold {i} overlaps"


def test_similarity_metric_is_the_indel_ratio():
    """Levenshtein.ratio normalises the INDEL distance: a substitution costs 2, not 1.

    Getting this wrong changes every similarity, every cluster and every split, so it is asserted
    rather than documented. 'abc' vs 'abd' needs one deletion plus one insertion: 1 - 2/6 = 2/3.
    """
    try:
        import Levenshtein
    except ImportError:
        return _skip("Levenshtein not installed")
    assert abs(Levenshtein.ratio("abc", "abd") - 2 / 3) < 1e-12, Levenshtein.ratio("abc", "abd")
    assert Levenshtein.distance("abc", "abd") == 1          # the edit distance itself is unchanged


def test_label_tensor_shape_and_sparsity():
    """174 x 75 x 84, 893,928 resolved pairs, 5,281 contacts (0.59 %) at 4.5 A."""
    if not DATA.exists():
        return _skip(f"{DATA.relative_to(ROOT)} is in the Zenodo archive; see DATA.md")
    z = np.load(DATA, allow_pickle=True)
    assert z["Y"].shape == (174, 75, 84), z["Y"].shape
    assert int(z["S"].sum()) == 893928, int(z["S"].sum())
    assert int(z["Y"].sum()) == 5281, int(z["Y"].sum())
    assert abs(100 * z["Y"].sum() / z["S"].sum() - 0.59) < 0.01
    assert ((z["Y"] == 1) & (z["S"] == 0)).sum() == 0, "a contact outside the resolved mask"
    rows = pd.Series([str(x) for x in z["row_imgt_region"]]).value_counts()
    assert dict(rows) == {"CDR3": 33, "CDR2": 15, "FR4": 14, "CDR1": 13}, dict(rows)
    cols = z["col_imgt_number"].astype(int)
    assert (cols <= 90).sum() == 41 and (cols >= 1050).sum() == 43


def test_stored_predictions_match_the_frozen_labels():
    """Every released prediction array is keyed to the same complexes, labels and masks."""
    if not PAIRMAPS.exists():
        return _skip(f"{PAIRMAPS.relative_to(ROOT)} is in the Zenodo archive; see DATA.md")
    md = pd.read_csv(ROOT / "data/processed/complex_metadata.tsv", sep="\t")
    lab = np.load(DATA, allow_pickle=True)
    for proto in ("random", "double", "allele", "temporal", "tcr70"):
        f = ROOT / f"results/tables/pairmaps_{proto}_p84v2.npz"
        if not f.exists():
            continue
        pm = np.load(f, allow_pickle=True)
        assert [str(x) for x in pm["pdb"]] == list(md.pdb), proto
        assert np.array_equal(pm["clu70"], md.clu70.values), proto
        assert np.array_equal(pm["Y"].astype(np.int8), lab["Y"]), proto
        assert np.array_equal(pm["S"].astype(np.int8), lab["S"]), proto
        # Complexes the protocol does not evaluate carry NaN by design; the coverage mask says
        # which rows are real, so range checks apply only there.
        cov = pm["covered"].astype(bool)
        for key in ("pred_model", "pred_freq"):
            p = pm[key]
            assert np.isfinite(p[cov]).all(), f"{proto}/{key}: non-finite inside coverage"
            assert np.isnan(p[~cov]).all(), f"{proto}/{key}: a value outside coverage"
            assert p[cov].min() >= 0 and p[cov].max() <= 1, f"{proto}/{key}: outside [0, 1]"
        assert int(cov.sum()) == {"random": 174, "double": 30, "allele": 64,
                                  "temporal": 60, "tcr70": 174}[proto], proto


def test_figure_1_is_the_approved_asset():
    """Figure 1 is the approved file, byte for byte, not a regenerated one."""
    import hashlib
    png = ROOT / "figures/main/Figure_1.png"
    if not png.exists():
        return _skip("figures/main/Figure_1.png not present")
    src = (ROOT / "figures/main/Figure_1.SOURCE.md").read_text()
    want = src.split("`figures/main/Figure_1.png` | `", 1)[1].split("`", 1)[0]
    got = hashlib.sha256(png.read_bytes()).hexdigest()
    assert got == want, f"Figure 1 is not the approved asset\n  want {want}\n  got  {got}"
    b = png.read_bytes()
    assert (int.from_bytes(b[16:20], "big"), int.from_bytes(b[20:24], "big")) == (6480, 4230)


def test_no_absolute_home_paths_in_shipped_code():
    """No advertised command may depend on the author's home directory."""
    # Built at runtime so this file does not itself contain the strings it looks for.
    needles = ["/" + "home" + "/", "C:" + chr(92) + "Users" + chr(92)]
    # figures/source/fig1_approved_20261102/ is the figure exactly as the author supplied it, kept
    # unmodified as the provenance record. The active source, figures/source/fig1_arial_20261102/,
    # is portable and is the one any command uses.
    archival = {"fig1_approved_20261102"}
    bad = []
    for d in ("scripts", "src", "tests", "md/analysis", "figures/source"):
        p = ROOT / d
        if not p.exists():
            continue
        for f in sorted(p.rglob("*.py")):
            if "__pycache__" in f.parts or archival & set(f.parts):
                continue
            t = f.read_text(errors="ignore")
            hit = next((n for n in needles if n in t), None)
            if hit:
                bad.append(f"{f.relative_to(ROOT)}: {hit!r}")
    assert not bad, "absolute home paths found:\n  " + "\n  ".join(bad)


def main() -> int:
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = skipped = 0
    for fn in fns:
        name = fn.__name__
        try:
            fn()
            print(f"[ ok ] {name}")
        except _Skip as e:
            skipped += 1
            print(f"[skip] {name}: {e}")
        except AssertionError as e:
            failed += 1
            print(f"[FAIL] {name}: {e}")
    passed = len(fns) - failed - skipped
    tail = f", {skipped} skipped (frozen data not present; see DATA.md)" if skipped else ""
    print(f"\n{passed}/{len(fns) - skipped} checks passed{tail}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
