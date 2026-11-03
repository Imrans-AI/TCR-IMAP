#!/usr/bin/env python3
"""p84v2 germline diagonal: executes scripts/40e_germline_diagonal_p82.py's own source (same
statistic, 2,000-resample bootstrap, seed 42, Holm across six loops) with its inputs pointed at the
p84v2 descriptive build and its output renamed. Validation first: the same exec with only the output
path redirected (to a temporary file) must reproduce results/tables/germline_diagonal_p82.csv exactly.
Output (refuses to overwrite): results/tables/germline_diagonal_p84v2.csv
usage: python scripts/53_germline_diagonal_p84v2.py
"""
import contextlib
import io
import sys
import tempfile
from pathlib import Path

import pandas as pd

SRC = Path("scripts/40e_germline_diagonal_p82.py").read_text()
OUT = "results/tables/germline_diagonal_p84v2.csv"
assert not Path(OUT).exists(), f"refusing to overwrite existing {OUT}"


def run(src, reps):
    for old, new, n in reps:
        assert src.count(old) == n, f"{old!r}: expected {n}, found {src.count(old)}"
        src = src.replace(old, new)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        exec(compile(src, "40e [p84v2 wrapper]", "exec"), {"__name__": "__main__"})
    return buf.getvalue()


tmp = Path(tempfile.mkdtemp()) / "germline_diagonal_p82_recheck.csv"
run(SRC, [("results/tables/germline_diagonal_p82.csv", str(tmp), 2)])
a, b = pd.read_csv(tmp), pd.read_csv("results/tables/germline_diagonal_p82.csv")
same = a.equals(b)
print(f"validation: 40e re-executed on the p82v1 inputs reproduces germline_diagonal_p82.csv exactly: {same}")
if not same:
    print("STOP"); sys.exit(1)
log = run(SRC, [("results/tables/contact_freq_gdomain180_descriptive.npz", "results/tables/contact_freq_gdomain182_descriptive_p84v2.npz", 1),
                ("results/tables/gdomain180_region_classification.csv", "results/tables/gdomain182_region_classification_p84v2.csv", 1),
                ("results/tables/germline_diagonal_p82.csv", OUT, 2)])
print(log)
