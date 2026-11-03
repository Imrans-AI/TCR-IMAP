#!/usr/bin/env python3
"""characterise.assembly_problem - the periodic-boundary guard (criteria document, 24 October 2026). Synthetic
chains with known answers; no trajectory is read.
usage: python md/analysis/tests/test_pbc_guard.py
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import characterise as ch   # noqa: E402

PASS, FAIL = [], []
def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name); print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' - ' + detail) if detail else ''}")

def chain(n, start, step=3.8):
    return np.array(start) + np.outer(np.arange(n), [step, 0, 0])

BOX = 150.0
A = chain(30, [0, 0, 0]); B = chain(30, [0, 6, 0]); C = chain(20, [0, -6, 0])     # three touching chains
X = np.vstack([A, B, C]); idx = {"A": np.arange(30), "B": np.arange(30, 60), "C": np.arange(60, 80)}
check("assembled complex passes", ch.assembly_problem(X, idx) is None)
Y = X.copy(); Y[idx["B"][12:]] += [0, 0, BOX]                                     # B broken across the box
r = ch.assembly_problem(Y, idx); check("a chain broken across the box is flagged", r is not None and "B broken" in r, str(r))
Z = X.copy(); Z[idx["C"]] += [BOX, 0, 0]                                          # C whole, in another image
r = ch.assembly_problem(Z, idx); check("a whole chain in another periodic image is flagged", r is not None and "not assembled" in r, str(r))
W = X.copy(); W[idx["A"][5]] += [0, 0, 1.2]                                       # a cis-like short/odd step, still bonded
check("ordinary geometry (step <= 4.5 A) passes", ch.assembly_problem(W, idx) is None)
print(f"\n{len(PASS)}/{len(PASS) + len(FAIL)} checks passed"); sys.exit(1 if FAIL else 0)
