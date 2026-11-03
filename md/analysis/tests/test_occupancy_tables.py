#!/usr/bin/env python3
"""occupancy_tables() on synthetic occupancy arrays with the real p84v2 label shapes, using the label
containers built by run_occupancy_analysis.main()'s OWN statements (extracted from its source and executed),
not containers this test builds for itself. Added 24 Oct 2026 after the authorised run crashed at
pipeline.py:51: every earlier test built its own containers as numpy arrays, so the run script's list
construction was never exercised. Reverting that construction to lists fails this test.

No trajectory and no score is read.
usage: python md/analysis/tests/test_occupancy_tables.py [path/to/run_occupancy_analysis.py]
       (the optional path lets the test be pointed at another version of the run script, e.g. the pre-fix one)
"""
import ast
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
A = HERE.parent
sys.path.insert(0, str(A))
from mdmap import ROOT, SYSTEMS          # noqa: E402
from pipeline import occupancy_tables    # noqa: E402

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' - ' + detail) if detail else ''}")


# ---- the run script's own construction of tl and mp, executed verbatim
RUN_SCRIPT = Path(sys.argv[1]) if len(sys.argv) > 1 else A / "run_occupancy_analysis.py"
print(f"  run script read: {RUN_SCRIPT}")
src = RUN_SCRIPT.read_text()
main = next(n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef) and n.name == "main")
stmts = [n for n in ast.walk(main) if isinstance(n, ast.Assign)
         and any(isinstance(t, ast.Tuple) and {getattr(e, "id", None) for e in t.elts} & {"tl", "mp"}
                 or getattr(t, "id", None) in ("tl", "mp") for t in n.targets)]
check("main() assigns tl and mp", {getattr(t, "id", None) for s in stmts for t in s.targets} >= {"tl", "mp"},
      "; ".join(ast.get_source_segment(src, s) for s in stmts))
lab = np.load(ROOT / "data/processed/pair_maps_p84v2.npz", allow_pickle=True)
ns = {"np": np, "lab": lab}
for s in stmts:
    exec(compile(ast.Module(body=[s], type_ignores=[]), "run_occupancy_analysis.py:main", "exec"), ns)
tl, mp = ns["tl"], ns["mp"]
print(f"  containers as main() builds them: tl {type(tl).__name__} {getattr(tl, 'dtype', '')} {len(tl)}, "
      f"mp {type(mp).__name__} {getattr(mp, 'dtype', '')} {len(mp)}")
check("tl and mp are numpy arrays (indexable by np.where)", isinstance(tl, np.ndarray) and isinstance(mp, np.ndarray))
check("their values equal the label file's row and column labels",
      list(map(str, tl)) == [str(x) for x in lab["row_label"]] and list(map(str, mp)) == [str(x) for x in lab["col_label"]])

# ---- synthetic Phase A results with the real label shapes (75 x 84), labels exactly as main() builds them
lpdb = [str(x) for x in lab["pdb"]]
labels = {p: (lab["Y"][lpdb.index(p)], lab["S"][lpdb.index(p)]) for p in SYSTEMS}
REPS = (1, 2, 3)
rng = np.random.default_rng(11)
occ = {}
for p in SYSTEMS:
    for r in REPS:
        o = rng.random(labels[p][0].shape)
        od = rng.random(labels[p][0].shape)
        odi = rng.random(labels[p][0].shape)
        occ[(p, r)] = dict(occ=o, occ_docked=od, occ_docked_identity=odi, n_window=9001,
                           n_docked=int(rng.integers(0, 9002)), n_docked_identity=int(rng.integers(0, 9002)),
                           ref_count=300, n_ref_pairs=310, persisting=rng.integers(0, 311, 9001),
                           counts=rng.integers(100, 400, 9001))
# one replicate with no docked occupancy at all, as run_trajectory returns when docking is off
occ[(SYSTEMS[0], 1)]["occ_docked"] = None

try:
    pairs, dock = occupancy_tables(occ, labels, tl, mp)
    check("occupancy_tables runs on the run script's containers", True)
except Exception as e:                                       # the 24 Oct crash lands here
    check("occupancy_tables runs on the run script's containers", False, f"{type(e).__name__}: {e}")
    print(f"\n{len(PASS)}/{len(PASS) + len(FAIL)} checks passed\nFAILED: " + ", ".join(FAIL))
    sys.exit(1)

n_res = {p: int(labels[p][1].sum()) for p in SYSTEMS}
check("one pair row per resolved pair per replicate", len(pairs) == 3 * sum(n_res.values()),
      f"{len(pairs):,} = 3 x {sum(n_res.values()):,}")
check("pair-table columns are the frozen schema",
      list(pairs.columns) == ["pdb", "rep", "tcr_slot", "mhc_position", "crystal_contact", "occupancy",
                              "occupancy_docked_count", "occupancy_docked_identity"], str(list(pairs.columns)))

# every row carries the right labels and values for its (row, column) cell
ok_lab = ok_val = ok_y = True
for (p, r), res in occ.items():
    Y, S = labels[p]
    ti, mj = np.where(S.astype(bool))
    sub = pairs[(pairs.pdb == p) & (pairs.rep == r)]
    ok_lab &= (sub.tcr_slot.values == tl[ti]).all() and (sub.mhc_position.values == mp[mj]).all()
    ok_y &= (sub.crystal_contact.values == Y[ti, mj]).all()
    ok_val &= np.array_equal(sub.occupancy.values, res["occ"][ti, mj]) \
        and np.array_equal(sub.occupancy_docked_identity.values, res["occ_docked_identity"][ti, mj]) \
        and (np.isnan(sub.occupancy_docked_count.values).all() if res["occ_docked"] is None
             else np.array_equal(sub.occupancy_docked_count.values, res["occ_docked"][ti, mj]))
check("tcr_slot / mhc_position are the labels of each row's own cell", ok_lab)
check("crystal_contact equals Y at each row's cell", ok_y)
check("occupancy columns equal the input arrays at each row's cell (NaN where docking is absent)", ok_val)
ins_rows = sorted(set(pairs.tcr_slot) - {x for x in pairs.tcr_slot if x[-1].isdigit()})
check("insertion-code TCR rows are labelled as such", set(ins_rows) <= {str(x) for x in tl}, ", ".join(ins_rows))
check("MHC 1061A and 1072A appear as their own columns", {"1061A", "1072A"} <= set(pairs.mhc_position))

check("one docked-frame row per replicate", len(dock) == len(occ))
for c in ("contributes_count", "contributes_identity", "readings_agree", "docked_fraction_count",
          "docked_fraction_identity", "persisting_median", "count_median"):
    check(f"docked-frame table has {c}", c in dock.columns)

print(f"\n{len(PASS)}/{len(PASS) + len(FAIL)} checks passed")
if FAIL:
    print("FAILED: " + ", ".join(FAIL))
sys.exit(1 if FAIL else 0)
