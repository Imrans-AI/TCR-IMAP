#!/usr/bin/env python3
"""End-to-end test of run_occupancy_analysis.main() ITSELF - its own label construction, its own phase
ordering, its own table writing and manifest - on the restrained NPT trajectories with SYNTHETIC scores.
Added 24 Oct 2026 after the authorised run crashed in main(): until then no test had executed main(), so a
defect in its own wiring could not be caught.

This is not the analysis. Inputs and outputs are redirected into a scratch tree:
  <scratch>/root/md/runs/<pdb>/rep<r>/prod.xtc  -> symlink to the real md/runs/<pdb>/rep<r>/npt.xtc (1 ns, 101 frames)
  <scratch>/root/md/runs/<pdb>/em/em.gro        -> symlink to the real em.gro (characterisation reference)
  <scratch>/root/results/tables/                -> every output table
  <scratch>/root/scores_synthetic.npz           -> random pred_model / pred_freq with the real pdb order, Y and S
The guards are replaced by a stub (the checkpoint check is meaningless for NPT), and the analysis window is set
to the NPT trajectory (0-1000 ps, 101 frames). Everything else in main() runs as written.

Hard safety, enforced by a Python audit hook that RAISES rather than logs:
  - no real production trajectory or checkpoint (md/runs/*/rep*/prod.*) is opened;
  - the real score file results/tables/pairmaps_tcr70_p84v2.npz is never opened;
  - nothing under the real results/tables/ is opened for writing.
Before and after: every real results/tables/md_* file (none before the authorised run; the run's outputs after it)
is byte-identical, and md/analysis/RUN_AUTHORISED is byte-identical.

A line tracer records which lines of the Phase A / A2 / B code executed, so coverage is measured rather
than assumed.

usage: python md/analysis/tests/run_main_on_npt.py <scratch_dir> [--min-docked N] [--run-script PATH]
  --min-docked N   lower stats.MIN_DOCKED_FRAMES for this test only, to exercise the contributing branch of
                   the docked-frame analysis (101 NPT frames never reach the production 1,000)
  --run-script P   execute main() from another version of the run script (e.g. the pre-fix one)
"""
import argparse
import hashlib
import importlib.util
import json
import os
import sys
import threading
import traceback
import warnings
from collections import defaultdict
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

A = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(A))
from mdmap import ROOT, SYSTEMS   # noqa: E402
import stats as st                # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("scratch")
ap.add_argument("--min-docked", type=int, default=None)
ap.add_argument("--run-script", default=str(A / "run_occupancy_analysis.py"))
args = ap.parse_args()

REAL_TABLES = (ROOT / "results/tables").resolve()
REAL_SCORES = (ROOT / "results/tables/pairmaps_tcr70_p84v2.npz").resolve()
REAL_RUNS = (ROOT / "md/runs").resolve()
AUTH = A / "RUN_AUTHORISED"
REPS = (1, 2, 3)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest() if Path(p).exists() else None


def real_md_outputs():
    """{name: sha256} of every real results/tables/md_* file. Before the authorised run this was empty; after it, the
    run's outputs exist. Either way the test must leave the set and every byte unchanged (checked after)."""
    return {x.name: sha(x) for x in sorted(REAL_TABLES.glob("md_*"))}


# ---------------- preconditions
pre_outputs, pre_auth = real_md_outputs(), sha(AUTH)
print(f"before: real results/tables md_* files: {len(pre_outputs)} (sha256 recorded); "
      f"RUN_AUTHORISED sha256 {str(pre_auth)[:16]}...")

# ---------------- the scratch tree
S = Path(args.scratch).resolve()
assert not str(S).startswith(str(ROOT.resolve())), "scratch must be outside the repository"
SR = S / "root"
if SR.exists():
    import shutil
    shutil.rmtree(SR)
for p in SYSTEMS:
    for r in REPS:
        src = ROOT / f"md/runs/{p}/rep{r}/npt.xtc"
        dst = SR / f"md/runs/{p}/rep{r}/prod.xtc"
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.symlink_to(src)
        assert dst.resolve().name == "npt.xtc", dst.resolve()
    em = SR / f"md/runs/{p}/em/em.gro"
    em.parent.mkdir(parents=True, exist_ok=True)
    em.symlink_to(ROOT / f"md/runs/{p}/em/em.gro")
OUT = SR / "results/tables"
OUT.mkdir(parents=True)
lab = np.load(ROOT / "data/processed/pair_maps_p84v2.npz", allow_pickle=True)
rng = np.random.default_rng(2024)
SCORES = SR / "scores_synthetic.npz"
np.savez(SCORES, pdb=lab["pdb"], Y=lab["Y"], S=lab["S"],
         pred_model=rng.random(lab["Y"].shape).astype(np.float32),
         pred_freq=rng.random(lab["Y"].shape).astype(np.float32))
print(f"scratch tree: {SR}  (12 prod.xtc -> npt.xtc symlinks, synthetic scores, empty output dir)")


# ---------------- the audit hook: raise on anything production
def _hook(event, a):
    if event != "open" or not a or not isinstance(a[0], (str, bytes, os.PathLike)):
        return
    path = os.path.realpath(os.fsdecode(a[0]))
    mode, flags = (a[1] if len(a) > 1 else None), (a[2] if len(a) > 2 else 0)
    writing = (isinstance(mode, str) and any(c in mode for c in "wax+")) or \
              (mode is None and isinstance(flags, int) and flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT))
    if path == str(REAL_SCORES):
        raise RuntimeError(f"SAFETY: attempted to open the real score file {path}")
    if path.startswith(str(REAL_RUNS)) and os.path.basename(path).startswith("prod."):
        raise RuntimeError(f"SAFETY: attempted to open a production MD file {path}")
    if writing and path.startswith(str(REAL_TABLES) + os.sep):
        raise RuntimeError(f"SAFETY: attempted to write under the real results/tables: {path}")


sys.addaudithook(_hook)

# ---------------- load the run script and redirect it
spec = importlib.util.spec_from_file_location("run_occupancy_analysis", args.run_script)
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)
R.HERE = A                     # residue maps, chain segments: the real ones
R.ROOT = SR                    # xtc_of, characterisation root (em.gro): the scratch tree
R.OUT = OUT
R.SCORES = SCORES
R.T_START_PS, R.T_END_PS, R.WINDOW_FRAMES = 0.0, 1000.0, 101
R.guards = lambda: {(p, r): R.NSTEPS for p in SYSTEMS for r in REPS}
assert R.LABELS.resolve() == (ROOT / "data/processed/pair_maps_p84v2.npz").resolve()
if args.min_docked is not None:
    st.MIN_DOCKED_FRAMES = args.min_docked
print(f"run script: {args.run_script}; window 0-1000 ps (101 frames); MIN_DOCKED_FRAMES {st.MIN_DOCKED_FRAMES}")

# ---------------- line tracer over the Phase A / A2 / B code
TARGETS = {str((A / f).resolve()): f for f in ("pipeline.py", "stats.py", "occupancy.py", "characterise.py",
                                               "contacts.py", "mdmap_p84v2.py")}
TARGETS[str(Path(args.run_script).resolve())] = "run_occupancy_analysis.py"
hits = defaultdict(set)


def _gtrace(frame, event, arg):
    fn = frame.f_code.co_filename
    if fn in TARGETS:
        hits[fn].add(frame.f_lineno)

        def _ltrace(frame, event, arg):
            if event == "line":
                hits[fn].add(frame.f_lineno)
            return _ltrace
        return _ltrace
    return None


sys.settrace(_gtrace)
threading.settrace(_gtrace)
crash = None
try:
    R.main()
except Exception as e:
    tb = traceback.extract_tb(e.__traceback__)[-1]
    crash = f"{type(e).__name__}: {e}  at {Path(tb.filename).name}:{tb.lineno}: {tb.line}"
finally:
    sys.settrace(None)
    threading.settrace(None)

# ---------------- postconditions, then the checks
post_outputs, post_auth = real_md_outputs(), sha(AUTH)
safety = post_outputs == pre_outputs and post_auth == pre_auth
print(f"\nafter: real results/tables md_* files byte-identical to before ({len(post_outputs)} files): "
      f"{post_outputs == pre_outputs}; RUN_AUTHORISED unchanged: {post_auth == pre_auth}")
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' - ' + detail) if detail else ''}")


check("safety: nothing written to the real results/tables, RUN_AUTHORISED untouched", safety)
check("main() completed without an exception", crash is None, crash or "")
if crash is None:
    written = sorted(x.name for x in OUT.glob("md_*"))
    check("all 12 output tables written (to scratch)", written == sorted(R.OUTPUTS), f"{len(written)} files")
    T = {n: pd.read_csv(OUT / n) for n in written}
    n_res = {p: int(lab["S"][[str(x) for x in lab["pdb"]].index(p)].sum()) for p in SYSTEMS}
    pairs = T["md_occupancy_pairs.csv"]
    check("pairs: one row per resolved pair per replicate", len(pairs) == 3 * sum(n_res.values()), f"{len(pairs):,}")
    # exactly the insertion-code labels that S resolves in these four complexes: MHC 1061A and 1072A in all four;
    # no CDR3 insertion row is resolved in any of them (MD_ANALYSIS_BUILD_REPORT.md), so none may appear
    lp = [str(x) for x in lab["pdb"]]
    want_m = {str(lab["col_label"][j]) for p in SYSTEMS for j in np.where(lab["S"][lp.index(p)].any(0))[0]
              if not str(lab["col_label"][j]).isdigit()}
    want_t = {str(lab["row_label"][i]) for p in SYSTEMS for i in np.where(lab["S"][lp.index(p)].any(1))[0]
              if str(lab["row_label"][i])[-1].isalpha()}
    got_m = {v for v in pairs.mhc_position.astype(str) if not v.isdigit()}
    got_t = {v for v in pairs.tcr_slot.astype(str) if v[-1].isalpha()}
    check("pairs: insertion-code labels are exactly those S resolves", got_m == want_m and got_t == want_t,
          f"MHC {sorted(got_m)}, TCR {sorted(got_t) or 'none'}")
    check("pairs: occupancy in [0, 1], no NaN", pairs.occupancy.between(0, 1).all() and pairs.occupancy.notna().all())
    man = dict(zip(T["md_occupancy_manifest.csv"]["item"], T["md_occupancy_manifest.csv"]["value"].astype(str)))
    check("manifest: pairs hash recorded before scores = the file written",
          man.get("md_occupancy_pairs_sha256_before_scores") == sha(OUT / "md_occupancy_pairs.csv"))
    check("manifest: characterisation hash recorded before scores = the file written",
          man.get("md_characterisation_summary_sha256_before_scores") == sha(OUT / "md_characterisation_summary.csv"))
    check("manifest: scores hash is the SYNTHETIC file", man.get("scores_sha256") == sha(SCORES))
    rho = T["md_occupancy_rho.csv"]
    check("rho: 34 rows = 2 scores x (4 systems x (1 + 3 reps) + pooled)", len(rho) == 34, str(len(rho)))
    check("rho: both scores, PRIMARY and CONTROL", set(rho.score) == {"TCR-IMAP", "pair-frequency prior"})
    check("rho: every value finite", rho.rho.notna().all())
    check("rho: synthetic scores give no association (|rho| < 0.1)", (rho.rho.abs() < 0.1).all(),
          f"max |rho| {rho.rho.abs().max():.3f}")
    strata = T["md_occupancy_strata.csv"]
    # a p-value is only meaningful where rho is defined; rows reported 'not tested' carry neither
    for name, tab in (("rho", rho), ("strata", strata), ("docked rho", T["md_occupancy_docked_rho.csv"])):
        bad = tab[tab.rho.isna() & tab.perm_p_approximate.notna()]
        check(f"{name}: no p-value reported for an undefined rho", bad.empty,
              "" if bad.empty else f"{len(bad)} rows with rho NaN but p = "
              + ", ".join(sorted({f'{v:.4f}' for v in bad.perm_p_approximate})))
    check("strata: 30 rows = 2 x (4 x 3 strata + 3 pooled)", len(strata) == 30, str(len(strata)))
    check("strata: the three fixed strata", set(strata.stratum) == {"4.5-6 A", "6-8 A", "> 8 A"})
    check("top-k: 8 rows", len(T["md_occupancy_topk.csv"]) == 8)
    check("check 1: 16 rows = 4 x (3 reps + system)", len(T["md_occupancy_check1_positive.csv"]) == 16)
    dr = T["md_occupancy_docked_rho.csv"]
    check("docked rho: 20 rows = 2 readings x 2 scores x (4 + pooled)", len(dr) == 20, str(len(dr)))
    check("docked rho: both readings", set(dr.reading) == {"count", "identity"})
    if st.MIN_DOCKED_FRAMES > 101:
        check("docked rho: no contributing replicate at the production threshold -> NaN, not an error",
              dr.rho.isna().all())
    else:
        sysrows = dr[dr.level == "system"]
        check("docked rho: contributing branch, every system has 3 replicates, values finite",
              (sysrows.n_replicates_contributing == 3).all() and dr.rho.notna().all())
    dock = T["md_occupancy_docked_frames.csv"]
    check("docked frames: 12 rows, both readings agree under restraint", len(dock) == 12 and dock.readings_agree.all())
    check("characterisation summary: 12 rows", len(T["md_characterisation_summary.csv"]) == 12)
    check("characterisation time series: 12 x (101 + 3 SASA) rows", len(T["md_characterisation_timeseries.csv"]) == 1248,
          str(len(T["md_characterisation_timeseries.csv"])))
    check("RMSF and H-bond tables non-empty", len(T["md_rmsf_per_residue.csv"]) > 0 and len(T["md_hbond_pairs.csv"]) > 0)

# ---------------- coverage, written for the report
cov = {TARGETS[f]: sorted(v) for f, v in hits.items()}
(S / f"coverage_{'min' + str(args.min_docked) if args.min_docked else 'prod'}.json").write_text(json.dumps(cov))
print(f"\nlines executed: " + ", ".join(f"{k} {len(v)}" for k, v in sorted(cov.items())))
print(f"\n{len(PASS)}/{len(PASS) + len(FAIL)} checks passed")
if FAIL:
    print("FAILED: " + ", ".join(FAIL))
sys.exit(1 if FAIL else 0)
