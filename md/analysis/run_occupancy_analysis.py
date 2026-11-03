#!/usr/bin/env python3
"""THE MD occupancy analysis - frozen; runs once, on production data, only when authorised.

Refuses to run unless
  1. md/analysis/RUN_AUTHORISED exists (created by the author, never by this code or by an assistant);
  2. all 12 production replicates are complete: checkpoint step of md/runs/<pdb>/rep<r>/prod.cpt is
     exactly 50,000,000 (read with `gmx dump -cp`, metadata only);
  3. no results/tables/md_occupancy_*.csv exists yet (the analysis is run once, never overwritten).

Then, in this order:
  Phase A - occupancy. For each system and replicate, one pass over prod.xtc: window 20-200 ns
            (9,001 frames, asserted), per-pair occupancy (minimum heavy-atom distance <= 4.5 A,
            minimum image), docked-frame occupancy and TCR-pMHC counts. Written to disk with hashes
            BEFORE any score is loaded.
  Phase B - statistics. Loads TCR-IMAP tcr70 out-of-fold pairwise scores (pred_model) and the
            within-fold pair-frequency prior (pred_freq) from results/tables/pairmaps_tcr70_p84v2.npz
            and computes every statistic in md/MD_SYSTEM_SELECTION_CRITERIA.md, "Analysis
            definitions - fixed 20 October 2026", plus check 1 and the 13 Oct docked-frame analysis.
Outputs (results/tables/):
  md_occupancy_pairs.csv          per pair x replicate occupancy (all + docked frames)       [phase A]
  md_occupancy_docked_frames.csv  per replicate docked fraction, contribution, counts        [phase A]
  md_occupancy_check1_positive.csv crystal-positive occupancy per replicate and system       [phase B]
  md_occupancy_rho.csv            Spearman rho per system (approx. permutation p), per
                                  replicate, pooled; TCR-IMAP and prior control              [phase B]
  md_occupancy_topk.csv           top-k crystal-negative fraction >= 0.25 vs base rate       [phase B]
  md_occupancy_docked_rho.csv     docked-frame secondary analysis                            [phase B]
  md_occupancy_manifest.csv       git commit, input hashes, versions, timings
usage (see results/review/MD_ANALYSIS_BUILD_REPORT.md for the exact command):
  python md/analysis/run_occupancy_analysis.py
"""
import hashlib
import platform
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
from mdmap import ROOT, SYSTEMS  # noqa: E402

REPS = (1, 2, 3)
NSTEPS = 50_000_000
T_START_PS, T_END_PS, WINDOW_FRAMES = 20_000.0, 200_000.0, 9_001
AUTH = HERE / "RUN_AUTHORISED"
OUT = ROOT / "results/tables"
OUTPUTS = ["md_occupancy_pairs.csv", "md_occupancy_docked_frames.csv", "md_occupancy_check1_positive.csv",
           "md_occupancy_rho.csv", "md_occupancy_topk.csv", "md_occupancy_docked_rho.csv",
           "md_occupancy_strata.csv", "md_characterisation_summary.csv",
           "md_characterisation_timeseries.csv", "md_rmsf_per_residue.csv", "md_hbond_pairs.csv",
           "md_occupancy_manifest.csv"]
# Characterisation settings (criteria document, 24 Oct 2026). SASA is subsampled because Shrake-Rupley
# costs orders of magnitude more per frame than the contact pass; the stride is justified by measurement
# in results/review/MD_DECISIONS_REPORT.md. Hydrogen bonds are computed on every window frame.
SASA_STRIDE = 50          # every 50th window frame = 1 ns spacing, 181 frames per replicate
HBOND_STRIDE = 1
# p84v2 artifacts: the corrected labels (174 x 75 x 84, insertion codes keyed) and the tcr70 out-of-fold
# predictions trained on them. The residue maps and chain segments are the p84v2 ones (gates G1-G3 re-passed
# 24 Oct 2026; see results/review/FINISH_PENDING_REPORT.md).
SCORES = ROOT / "results/tables/pairmaps_tcr70_p84v2.npz"
LABELS = ROOT / "data/processed/pair_maps_p84v2.npz"
GMXRC = "/usr/local/gromacs-2023.2/bin/GMXRC"


def checkpoint_step(cpt):
    """Last completed step recorded in a checkpoint (metadata only), as md/lib.sh laststep()."""
    cmd = (f"command -v gmx >/dev/null 2>&1 || source {GMXRC}; "
           f"gmx dump -cp '{cpt}' 2>/dev/null | awk '/^ *step *=/{{print $3; exit}}'")
    out = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True).stdout.strip()
    return int(out) if out.isdigit() else -1


def refuse(msg):
    print(f"REFUSED: {msg}")
    sys.exit(2)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def guards():
    if not AUTH.exists():
        refuse(f"{AUTH.relative_to(ROOT)} does not exist; the author creates it to authorise the single run")
    status = {(p, r): checkpoint_step(ROOT / f"md/runs/{p}/rep{r}/prod.cpt") for p in SYSTEMS for r in REPS}
    bad = {k: v for k, v in status.items() if v != NSTEPS}
    if bad:
        refuse("production incomplete: " + ", ".join(f"{p} rep{r} step {s}" for (p, r), s in bad.items()))
    exist = [o for o in OUTPUTS if (OUT / o).exists()]
    if exist:
        refuse("outputs already exist (the analysis runs once): " + ", ".join(exist))
    return status


def main():
    t0 = time.time()
    status = guards()
    print("guards passed: RUN_AUTHORISED present; all 12 checkpoints at step 50,000,000; no prior outputs")

    from pipeline import (occupancy_phase, occupancy_tables, statistics_phase, characterisation_phase,
                          PRIMARY, CONTROL)

    from mdmap_p84v2 import groups_from_map_p84v2

    lab = np.load(LABELS, allow_pickle=True)
    lpdb = [str(x) for x in lab["pdb"]]
    # numpy arrays, not lists: occupancy_tables indexes them with np.where arrays (crash fixed 24 Oct 2026)
    tl = np.asarray(lab["row_label"], dtype=str); mp = np.asarray(lab["col_label"], dtype=str)
    labels = {p: (lab["Y"][lpdb.index(p)], lab["S"][lpdb.index(p)]) for p in SYSTEMS}
    rmaps = {p: pd.read_csv(HERE / f"residue_map_{p}_p84v2.tsv", sep="\t") for p in SYSTEMS}
    segs = pd.read_csv(HERE / "chain_segments_p84v2.tsv", sep="\t")
    xtc_of = lambda p, r: ROOT / f"md/runs/{p}/rep{r}/prod.xtc"

    # ---------------- Phase A: occupancy, no scores in memory ----------------
    print("phase A: occupancy (20-200 ns, 9,001 frames per replicate)")
    tA = time.time()
    occ = occupancy_phase(SYSTEMS, REPS, xtc_of, rmaps, segs, T_START_PS, T_END_PS, WINDOW_FRAMES,
                          groups_fn=groups_from_map_p84v2)
    pairs, dock = occupancy_tables(occ, labels, tl, mp)
    pairs.to_csv(OUT / "md_occupancy_pairs.csv", index=False)
    dock.to_csv(OUT / "md_occupancy_docked_frames.csv", index=False)
    phaseA_hash = sha256(OUT / "md_occupancy_pairs.csv")
    tA = time.time() - tA
    print(f"phase A written ({tA / 60:.1f} min); md_occupancy_pairs.csv sha256 {phaseA_hash}")

    # ---------------- Phase A2: characterisation, still score-free ----------------
    print("phase A2: simulation quality and interface characterisation (no scores)")
    tA2 = time.time()
    summ, series, rmsf, hbp = characterisation_phase(SYSTEMS, REPS, xtc_of, rmaps, segs,
                                                     T_START_PS, T_END_PS, ROOT,
                                                     sasa_stride=SASA_STRIDE, hbond_stride=HBOND_STRIDE)
    summ.to_csv(OUT / "md_characterisation_summary.csv", index=False)
    series.to_csv(OUT / "md_characterisation_timeseries.csv", index=False)
    rmsf.to_csv(OUT / "md_rmsf_per_residue.csv", index=False)
    hbp.to_csv(OUT / "md_hbond_pairs.csv", index=False)
    charA_hash = sha256(OUT / "md_characterisation_summary.csv")
    tA2 = time.time() - tA2
    print(f"phase A2 written ({tA2 / 60:.1f} min); md_characterisation_summary.csv sha256 {charA_hash}")

    # crystal-frame distances for the stratified rho - structural, no score involved
    from mdmap_p84v2 import crystal_distance_matrix
    dists = {p: crystal_distance_matrix(p) for p in SYSTEMS}

    # ---------------- Phase B: scores loaded only now ----------------
    pm = np.load(SCORES, allow_pickle=True)
    assert list(pm["pdb"]) == lpdb, "score/label complex order differs"
    scores = {}
    for name, key in ((PRIMARY, "pred_model"), (CONTROL, "pred_freq")):
        scores[name] = {}
        for p in SYSTEMS:
            i = lpdb.index(p)
            assert (pm["Y"][i] == lab["Y"][i]).all() and (pm["S"][i] == lab["S"][i]).all(), f"{p}: Y/S differ"
            scores[name][p] = pm[key][i].astype(float)
    tB = time.time()
    res = statistics_phase(occ, labels, scores, SYSTEMS, REPS, dists=dists)
    tB = time.time() - tB
    res["check1"].to_csv(OUT / "md_occupancy_check1_positive.csv", index=False)
    res["rho"].to_csv(OUT / "md_occupancy_rho.csv", index=False)
    res["topk"].to_csv(OUT / "md_occupancy_topk.csv", index=False)
    res["docked_rho"].to_csv(OUT / "md_occupancy_docked_rho.csv", index=False)
    res["strata"].to_csv(OUT / "md_occupancy_strata.csv", index=False)

    git = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "md/analysis"],
                           capture_output=True, text=True).stdout.strip()
    import MDAnalysis, scipy
    man = [dict(item="git_commit", value=git), dict(item="md_analysis_tree_clean", value=not dirty),
           dict(item="python", value=platform.python_version()), dict(item="MDAnalysis", value=MDAnalysis.__version__),
           dict(item="numpy", value=np.__version__), dict(item="scipy", value=scipy.__version__),
           dict(item="labels_sha256", value=sha256(LABELS)), dict(item="scores_sha256", value=sha256(SCORES)),
           dict(item="md_occupancy_pairs_sha256_before_scores", value=phaseA_hash),
           dict(item="permutation", value="within-system, 10,000, seed 42, two-sided, approximate"),
           dict(item="docked_readings", value="count-based (primary) and identity-based, both reported"),
           dict(item="distance_strata", value="crystal-frame min heavy-atom distance: 4.5-6, 6-8, > 8 A"),
           dict(item="sasa", value=f"Shrake-Rupley, probe 1.4 A, 256 sphere points, stride {SASA_STRIDE}"),
           dict(item="hbond", value=f"d(D-A) <= 3.5 A, angle(D-H...A) >= 150 deg, stride {HBOND_STRIDE}"),
           dict(item="md_characterisation_summary_sha256_before_scores", value=charA_hash),
           dict(item="phase_A_minutes", value=round(tA / 60, 2)),
           dict(item="phase_A2_minutes", value=round(tA2 / 60, 2)),
           dict(item="phase_B_minutes", value=round(tB / 60, 2)),
           dict(item="total_minutes", value=round((time.time() - t0) / 60, 2))]
    for p in SYSTEMS:
        man.append(dict(item=f"residue_map_{p}_p84v2_sha256", value=sha256(HERE / f"residue_map_{p}_p84v2.tsv")))
        for r in REPS:
            x = xtc_of(p, r)
            man.append(dict(item=f"{p}_rep{r}_checkpoint_step", value=status[(p, r)]))
            man.append(dict(item=f"{p}_rep{r}_prod_xtc_bytes", value=x.stat().st_size))
    pd.DataFrame(man).to_csv(OUT / "md_occupancy_manifest.csv", index=False)

    with pd.option_context("display.width", 200, "display.max_columns", 20):
        for k in ("check1", "rho", "topk", "docked_rho", "strata"):
            print(f"\n--- {k} ---\n{res[k].round(4).to_string(index=False)}")
        print(f"\n--- docked frames ---\n{dock.round(4).to_string(index=False)}")
    print(f"\ndone in {(time.time() - t0) / 60:.1f} min; tables in results/tables/md_occupancy_*.csv")


if __name__ == "__main__":
    main()
