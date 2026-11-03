#!/usr/bin/env python3
"""Recompute ONLY the characterisation phase with periodic-boundary handling (criteria document, "Characterisation
recomputed with periodic-boundary handling", 24 October 2026).

Score-free: reads no score file, computes no occupancy and no correlation. The eight occupancy tables and the run
manifest are hashed before and after and must be byte-identical; the script stops otherwise.

  coordinate measures (backbone RMSD, interface RMSD, per-loop RMSF, buried SASA): the whole, jump-free
      trajectories written by md/analysis/pbc.py (qc_replicate.sh's procedure, every frame)
  hydrogen bonds: the original prod.xtc (periodic-aware), exactly as in the run, so they must not move

Window, SASA stride, H-bond stride, selections and references are imported from run_occupancy_analysis.py, so they
are the run's own. The superseded tables are copied to results/review/characterisation_fix/backup_superseded/ first;
if that copy already exists it must hold the run's original summary (hash in the run manifest) and is not overwritten.
Writes the four characterisation tables plus results/tables/md_characterisation_manifest.csv.

usage: python md/analysis/recompute_characterisation.py <dir with <pdb>_rep<r>/nj.xtc>
"""
import hashlib
import platform
import shutil
import subprocess
import sys
import time
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from mdmap import ROOT, SYSTEMS                                         # noqa: E402
import run_occupancy_analysis as RUN                                    # constants only; main() is not called
from pipeline import characterisation_phase                             # noqa: E402
from pbc import xtc_frames                                              # noqa: E402

T = ROOT / "results/tables"
OCC = ["md_occupancy_pairs.csv", "md_occupancy_rho.csv", "md_occupancy_strata.csv", "md_occupancy_topk.csv",
       "md_occupancy_check1_positive.csv", "md_occupancy_docked_frames.csv", "md_occupancy_docked_rho.csv",
       "md_occupancy_manifest.csv"]
CHAR = ["md_characterisation_summary.csv", "md_characterisation_timeseries.csv", "md_rmsf_per_residue.csv",
        "md_hbond_pairs.csv"]
BK = ROOT / "results/review/characterisation_fix/backup_superseded"
REPS = RUN.REPS


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    WD = Path(sys.argv[1]).resolve()
    t0 = time.time()
    occ_before = {n: sha(T / n) for n in OCC}
    man = pd.read_csv(T / "md_occupancy_manifest.csv").set_index("item")["value"]
    orig_summary = man["md_characterisation_summary_sha256_before_scores"]

    # ---- backups of the superseded tables (never overwritten once made)
    if BK.exists():
        assert sha(BK / "md_characterisation_summary.csv") == orig_summary, "backup is not the run's original"
        print(f"backups already present and verified original: {BK.relative_to(ROOT)}")
    else:
        assert sha(T / "md_characterisation_summary.csv") == orig_summary, \
            "results/tables summary is not the run's original - refusing to back up something else"
        BK.mkdir(parents=True)
        for n in CHAR:
            shutil.copy2(T / n, BK / n)
            assert sha(BK / n) == sha(T / n)
        print(f"superseded tables backed up to {BK.relative_to(ROOT)}")
    superseded = {n: sha(BK / n) for n in CHAR}

    # ---- inputs
    nj = lambda p, r: WD / f"{p}_rep{r}/nj.xtc"
    prod = lambda p, r: ROOT / f"md/runs/{p}/rep{r}/prod.xtc"
    frames = {}
    for p in SYSTEMS:
        for r in REPS:
            n, first, last = xtc_frames(nj(p, r), WD)
            assert n == 10001 and float(last) == 200000.0, f"{p} rep{r}: {n} frames, last {last}"
            frames[(p, r)] = n
    rmaps = {p: pd.read_csv(HERE / f"residue_map_{p}_p84v2.tsv", sep="\t") for p in SYSTEMS}
    segs = pd.read_csv(HERE / "chain_segments_p84v2.tsv", sep="\t")
    print(f"inputs: 12 whole + nojump trajectories (10,001 frames each); H-bonds from the original prod.xtc; "
          f"window {RUN.T_START_PS:.0f}-{RUN.T_END_PS:.0f} ps, SASA stride {RUN.SASA_STRIDE}, H-bond stride {RUN.HBOND_STRIDE}")

    # ---- the characterisation phase, nothing else
    summ, series, rmsf, hbp = characterisation_phase(SYSTEMS, REPS, nj, rmaps, segs, RUN.T_START_PS, RUN.T_END_PS,
                                                     ROOT, sasa_stride=RUN.SASA_STRIDE, hbond_stride=RUN.HBOND_STRIDE,
                                                     hbond_xtc_of=prod)
    summ.to_csv(T / "md_characterisation_summary.csv", index=False)
    series.to_csv(T / "md_characterisation_timeseries.csv", index=False)
    rmsf.to_csv(T / "md_rmsf_per_residue.csv", index=False)
    hbp.to_csv(T / "md_hbond_pairs.csv", index=False)

    occ_after = {n: sha(T / n) for n in OCC}
    unchanged = occ_before == occ_after
    for n in OCC:
        print(f"  {'unchanged' if occ_before[n] == occ_after[n] else 'CHANGED  '}  {n}  {occ_after[n][:16]}")
    assert unchanged, "an occupancy table changed - STOP"

    import MDAnalysis
    git = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "md/analysis"],
                           capture_output=True, text=True).stdout.strip()
    gmxv = subprocess.run(["bash", "-c", "command -v gmx >/dev/null || source /usr/local/gromacs-2023.2/bin/GMXRC; "
                                         "gmx --version 2>/dev/null | grep -m1 'GROMACS version'"],
                          capture_output=True, text=True).stdout.strip()
    rows = [dict(item="git_commit", value=git), dict(item="md_analysis_tree_clean", value=not dirty),
            dict(item="reason", value="periodic-boundary handling (criteria document, 24 October 2026)"),
            dict(item="coordinate_measures_trajectory", value="gmx convert-tpr SOLU; trjconv -pbc whole; trjconv -pbc nojump; every frame (md/analysis/pbc.py)"),
            dict(item="hbond_trajectory", value="original prod.xtc (periodic-aware), as in the run"),
            dict(item="assembly_guard", value="every window frame: consecutive Ca <= 4.5 A per chain, chains connected at 15 A"),
            dict(item="gromacs", value=gmxv), dict(item="python", value=platform.python_version()),
            dict(item="MDAnalysis", value=MDAnalysis.__version__),
            dict(item="window_ps", value=f"{RUN.T_START_PS:.0f}-{RUN.T_END_PS:.0f}"),
            dict(item="sasa_stride", value=RUN.SASA_STRIDE), dict(item="hbond_stride", value=RUN.HBOND_STRIDE),
            dict(item="occupancy_tables_unchanged", value=unchanged),
            dict(item="minutes", value=round((time.time() - t0) / 60, 2))]
    for n in CHAR:
        rows.append(dict(item=f"{n}_sha256", value=sha(T / n)))
        rows.append(dict(item=f"{n}_superseded_sha256", value=superseded[n]))
    for n in OCC:
        rows.append(dict(item=f"{n}_sha256_unchanged", value=occ_after[n]))
    for p in SYSTEMS:
        for r in REPS:
            rows.append(dict(item=f"{p}_rep{r}_prod_xtc_bytes", value=prod(p, r).stat().st_size))
            rows.append(dict(item=f"{p}_rep{r}_whole_nojump_frames", value=frames[(p, r)]))
    pd.DataFrame(rows).to_csv(T / "md_characterisation_manifest.csv", index=False)
    print(f"done in {(time.time() - t0) / 60:.1f} min; four characterisation tables rewritten, manifest written; "
          f"eight occupancy tables byte-identical")


if __name__ == "__main__":
    main()
