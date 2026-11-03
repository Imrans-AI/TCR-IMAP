#!/usr/bin/env python3
"""Smoke test of the full pipeline wiring WITHOUT production data and WITHOUT real scores:
pipeline.occupancy_phase + statistics_phase on the 12 restrained NPT trajectories (rep1-3, 1 ns,
101 frames each) with SYNTHETIC random scores, writing only to a scratch directory. The run script
itself is not executed; only its checkpoint reader is imported and called on checkpoint metadata.
usage: python md/analysis/tests/smoke_pipeline_npt.py <scratch_dir>
"""
import sys
import tempfile
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from mdmap import ROOT, SYSTEMS  # noqa: E402
import stats as st  # noqa: E402
from pipeline import occupancy_phase, occupancy_tables, statistics_phase  # noqa: E402
from run_occupancy_analysis import checkpoint_step, NSTEPS  # noqa: E402  (helper only; main() not called)

out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(tempfile.mkdtemp())
out.mkdir(parents=True, exist_ok=True)
REPS = (1, 2, 3)
lab = np.load(ROOT / "data/processed/pair_maps_p82v1.npz", allow_pickle=True)
lpdb = list(lab["pdb"])
labels = {p: (lab["Y"][lpdb.index(p)], lab["S"][lpdb.index(p)]) for p in SYSTEMS}
rmaps = {p: pd.read_csv(HERE / f"residue_map_{p}.tsv", sep="\t") for p in SYSTEMS}
segs = pd.read_csv(HERE / "chain_segments.tsv", sep="\t")

occ = occupancy_phase(SYSTEMS, REPS, lambda p, r: ROOT / f"md/runs/{p}/rep{r}/npt.xtc", rmaps, segs,
                      0.0, 1000.0, 101)
pairs, dock = occupancy_tables(occ, labels, lab["tcr_position_labels"], lab["mhc_positions"])
assert len(pairs) == sum(int(labels[p][1].sum()) for p in SYSTEMS) * 3
assert pairs.occupancy.between(0, 1).all() and not pairs.occupancy.isna().any()
# 101 frames < 1,000: none contribute. Column renamed 24 Oct 2026 (bfc76b6, decision 3: both docked readings);
# the same expectation now applies to both readings.
assert (~dock.contributes_count).all() and (~dock.contributes_identity).all()

rng = np.random.default_rng(7)
synth = {"synthetic A": {p: rng.random((84, 82)) for p in SYSTEMS},
         "synthetic B": {p: rng.random((84, 82)) for p in SYSTEMS}}
res = statistics_phase(occ, labels, synth, SYSTEMS, REPS, n_perm=999)
rho = res["rho"]
assert len(rho) == 2 * (4 * (1 + 3) + 1) and rho.rho.notna().all()
assert (res["topk"].k.values == np.tile([int((labels[p][0] & labels[p][1]).sum()) for p in SYSTEMS], 2)).all()
assert res["docked_rho"].rho.isna().all()                      # no contributing replicate -> NaN, not an error
assert len(res["check1"]) == 4 * 4

st.MIN_DOCKED_FRAMES = 50                                      # exercise the contributing branch only here
res2 = statistics_phase(occ, labels, synth, SYSTEMS, REPS, n_perm=999)
assert res2["docked_rho"].rho.notna().all() and (res2["docked_rho"].query("level == 'system'")
                                                  .n_replicates_contributing == 3).all()
st.MIN_DOCKED_FRAMES = 1000

for k, v in res.items():
    v.to_csv(out / f"smoke_{k}.csv", index=False)
pairs.to_csv(out / "smoke_pairs.csv", index=False); dock.to_csv(out / "smoke_docked_frames.csv", index=False)
with pd.option_context("display.width", 200):
    print(res["check1"].round(3).to_string(index=False))
    print(rho[rho.level != "replicate"].round(3).to_string(index=False))
    print(res["topk"].round(3).to_string(index=False))

# checkpoint reader (metadata only) on every production checkpoint, and on a missing file
st_rows = {(p, r): checkpoint_step(ROOT / f"md/runs/{p}/rep{r}/prod.cpt") for p in SYSTEMS for r in REPS}
assert checkpoint_step(ROOT / "md/runs/none/prod.cpt") == -1
done = sum(v == NSTEPS for v in st_rows.values())
print("\nproduction checkpoints: " + ", ".join(f"{p} rep{r} {v:,}" for (p, r), v in st_rows.items()))
print(f"complete: {done}/12 -> the run script would {'PASS' if done == 12 else 'REFUSE'} the completion guard")
print(f"\nsmoke test PASS; outputs in {out}")
