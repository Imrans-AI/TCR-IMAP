#!/usr/bin/env python3
"""Derive Results Tables 1 and 2 of the v4 manuscript from the committed artifacts.

Every numeric cell is read from a frozen file; none is typed by hand. The qualitative columns -
what a protocol controls and what role it plays - are editorial text and are written here so that
they sit beside the numbers they describe rather than in the builder.

Writes docs/manuscript/results_tables_v4.json, which build_manuscript_docx_v4.js reads. Re-run it
after any change to the underlying tables; the builder asserts the file is present and that its
provenance hashes still match.

Sources:
  data/splits_v2/protocols.json                      folds and complexes evaluated
  data/processed/complex_metadata.tsv                cohort size
  results/review/figures_p84v2/stats_20k_seed42.tsv  differences, intervals and p values

usage: python scripts/71_results_tables_v4.py
"""
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs/manuscript/results_tables_v4.json"

PROTO = ROOT / "data/splits_v2/protocols.json"
META = ROOT / "data/processed/complex_metadata.tsv"
STATS = ROOT / "results/review/figures_p84v2/stats_20k_seed42.tsv"
MIN_FOLD = 3                      # folds below this are dropped at evaluation
MINUS = "−"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def signed(x: float, nd: int = 4) -> str:
    return f"{x:+.{nd}f}".replace("-", MINUS)


def ci(lo: float, hi: float) -> str:
    return f"[{signed(lo)}, {signed(hi)}]"


def pstr(p: float) -> str:
    """The bootstrap floor is 1/20,000; a value at the floor is reported as p < 0.001."""
    return "< 0.001" if p <= 0.0005 else f"{p:.3f}"


# ------------------------------------------------------------------ Table 1
proto = json.loads(PROTO.read_text())
md = pd.read_csv(META, sep="\t")
N = len(md)
assert N == 174, N

# Display name, frozen key, and the two control columns. The qualifications here are the ones
# Methods 2.3 carries in full; they are repeated in short form so the table is readable alone.
T1_SPEC = [
    ("random", "random",
     "none; folds are drawn over complexes",
     "none",
     "reference point; the inflation other protocols remove is measured against it"),
    ("intersection (legacy)", "double",
     "partial; test complexes are selected at the intersection of a held-out cluster and a "
     "held-out allele group, but all non-test complexes remain in training, so the held-out groups "
     "are still seen",
     "partial, in the same sense",
     "legacy, retained for continuity with earlier work; not a double holdout and not a point on "
     "any strictness ranking"),
    ("allele-annotation holdout", "allele",
     "none",
     "holdout is by parsed allele annotation, not by MHC sequence: 108 complexes carry no "
     "annotation and never leave training, and for 32 of the 64 evaluated complexes an identical "
     "84-residue MHC window is present in training",
     "MHC-side redundancy, as far as the available annotation allows"),
    ("temporal", "temporal",
     "none",
     "none",
     "prospective-style split by deposition date"),
    ("tcr70", "tcr70",
     "whole TCR clusters at a similarity ratio of 0.70 are held out; zero train/test cluster "
     "overlap in every fold",
     "none; it partitions on TCR sequence only",
     "primary protocol; the predeclared comparison is made under it"),
]

table1 = []
for name, key, tcr_ctl, mhc_ctl, role in T1_SPEC:
    folds = proto[key]
    evaluated = sum(len(f) for f in folds if len(f) >= MIN_FOLD)
    dropped = sum(1 for f in folds if len(f) < MIN_FOLD)
    table1.append({
        "protocol": name,
        "folds": str(len(folds)) + (f" ({dropped} dropped, n < {MIN_FOLD})" if dropped else ""),
        "evaluated": f"{evaluated} of {N}",
        "tcr_control": tcr_ctl,
        "mhc_control": mhc_ctl,
        "role": role,
    })

# ------------------------------------------------------------------ Table 2
stats = pd.read_csv(STATS, sep="\t", float_precision="round_trip")


def row(group: str, comparison: str):
    r = stats[(stats.group == group) & (stats.comparison == comparison)]
    assert len(r) == 1, (group, comparison, len(r))
    return r.iloc[0]


T2_SPEC = [
    ("per-residue", "tcr70: TCR-IMAP vs Struct-TCS", "TCR-IMAP vs Struct-TCS",
     "per-residue", "tcr70",
     "**predeclared primary comparison.** Unresolved: the interval spans a disadvantage and an "
     "advantage more than twice the estimated difference"),
    ("pairwise", "tcr70: TCR-IMAP vs pair-frequency prior", "TCR-IMAP vs pair-frequency prior",
     "pairwise", "tcr70",
     "secondary, unadjusted; TCR-IMAP higher on 112 of 174 complexes"),
    ("pairwise", "tcr70: TCR-IMAP vs TEIM-Res", "TCR-IMAP vs TEIM-Res",
     "pairwise", "tcr70",
     "secondary, unadjusted; TCR-IMAP higher on 143 of 174 complexes"),
    ("per-residue", "random: TCR-IMAP vs Struct-TCS", "TCR-IMAP vs Struct-TCS",
     "per-residue", "random",
     "secondary, unadjusted; no difference is detected where train/test similarity is "
     "unrestricted"),
]

table2 = []
for group, comp_key, label, task, protocol, interp in T2_SPEC:
    r = row(group, comp_key)
    table2.append({
        "comparison": label,
        "task": task,
        "protocol": protocol,
        "delta": signed(float(r["diff"])),
        "ci": ci(float(r["lo"]), float(r["hi"])),
        "p": pstr(float(r["p"])),
        "n": int(r["n"]),
        "interpretation": interp,
    })

# the win counts quoted in the interpretation column must come from the same file
for t2, key in zip(table2[1:3], ("tcr70: TCR-IMAP vs pair-frequency prior",
                                 "tcr70: TCR-IMAP vs TEIM-Res")):
    wins = str(row("pairwise", key)["wins"])          # recorded as "112/174"
    w, tot = wins.split("/")
    assert f"higher on {w} of {tot}" in t2["interpretation"], (key, wins)

payload = {
    "generated_by": "scripts/71_results_tables_v4.py",
    "sources": {str(p.relative_to(ROOT)): sha(p) for p in (PROTO, META, STATS)},
    "table1": table1,
    "table2": table2,
}
OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")

print(f"wrote {OUT.relative_to(ROOT)}")
print("\nTable 1")
for r in table1:
    print(f"  {r['protocol']:<26} folds {r['folds']:<22} evaluated {r['evaluated']}")
print("\nTable 2")
for r in table2:
    print(f"  {r['comparison']:<34} {r['task']:<12} {r['protocol']:<8} "
          f"{r['delta']} {r['ci']}  p {r['p']}  n {r['n']}")
