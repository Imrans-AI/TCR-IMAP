#!/usr/bin/env python3
"""One authoritative uncertainty calculation for the identity-threshold sweep.

Diagnosed in results/review/SWEEP_UNCERTAINTY_RECONCILIATION.md: the sweep's committed CIs and p
values were produced by cluster_bootstrap_diff's DEFAULTS (2,000 resamples, seed 0) because the
sweep stored no per-complex arrays, while every other analysis in the project uses 20,000
resamples at seed 42. The estimand is the same and the point estimates are identical; the
difference is Monte-Carlo noise, and at 2,000 resamples the 0.70 interval's lower endpoint is
unstable around zero.

Per-complex arrays for all seven thresholds were recovered by scripts/62 under a validation gate
that reproduced the frozen aggregates exactly. This script recomputes every threshold at the
project standard, computes p directly (never inferred from whether the CI covers zero), and
writes an audit table carrying the superseded 2,000-resample values beside the authoritative ones.

Also, now that threshold 1.01 per-complex values exist, the size effect (unrestricted versus
size-matched random training at the 0.70-matched size) is given a proper PAIRED interval instead
of a difference of point estimates.

Writes results/tables/similarity_curve_p84v2_authoritative.csv and
results/review/sweep_authoritative/. The frozen similarity_curve_p84v2.csv is NOT modified.

usage: python scripts/69_sweep_authoritative_stats.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, ".")
from src.eval.metrics import cluster_bootstrap_diff

ROOT = Path(__file__).resolve().parent.parent
FR = ROOT/"results/review/filtered_rerun/filtered_per_complex.tsv"
SM = ROOT/"results/review/size_matched_control/per_complex_metrics.tsv"
FROZEN = ROOT/"results/tables/similarity_curve_p84v2.csv"
OUT = ROOT/"results/review/sweep_authoritative"; OUT.mkdir(parents=True, exist_ok=True)
STD = dict(n_boot=20000, seed=42)
OLD = dict(n_boot=2000, seed=0)          # the historical defaults


def boot(a, b, clu, **kw):
    d, lo, hi, p = cluster_bootstrap_diff(a, b, clu, **kw)
    return float(d), float(lo), float(hi), float(p)


def main():
    fr = pd.read_csv(FR, sep="\t")
    frozen = pd.read_csv(FROZEN).set_index("threshold")
    rows = []
    for thr in sorted(fr.threshold.unique(), reverse=True):
        g = fr[fr.threshold == thr]
        a, b, c = g.model_auprc.values, g.freq_auprc.values, g.clu70.values
        d_new, lo_n, hi_n, p_n = boot(a, b, c, **STD)
        d_old, lo_o, hi_o, p_o = boot(a, b, c, **OLD)
        fz = frozen.loc[thr]
        rows.append(dict(
            threshold=thr, n=len(g),
            model=round(float(np.nanmean(a)), 4), freq=round(float(np.nanmean(b)), 4),
            diff=round(d_new, 4),
            lo=round(lo_n, 4), hi=round(hi_n, 4), p=round(p_n, 4),
            ci_excludes_zero=bool(lo_n > 0 or hi_n < 0),
            superseded_lo=round(lo_o, 4), superseded_hi=round(hi_o, 4), superseded_p=round(p_o, 4),
            superseded_ci_excludes_zero=bool(lo_o > 0 or hi_o < 0),
            frozen_csv_lo=fz.lo, frozen_csv_hi=fz.hi, frozen_csv_p=fz.p,
            reproduces_frozen=bool(round(lo_o, 4) == round(float(fz.lo), 4)
                                   and round(hi_o, 4) == round(float(fz.hi), 4)
                                   and round(p_o, 4) == round(float(fz.p), 4))))
    t = pd.DataFrame(rows).sort_values("threshold", ascending=False)
    t.to_csv(OUT/"sweep_authoritative.tsv", sep="\t", index=False)
    # a drop-in replacement for the figure, same columns as the frozen CSV
    auth = t[["threshold", "n", "model", "freq", "diff", "lo", "hi", "p"]].copy()
    auth["mean_max_sim"] = [frozen.loc[x].mean_max_sim for x in auth.threshold]
    auth.sort_values("threshold").to_csv(ROOT/"results/tables/similarity_curve_p84v2_authoritative.csv", index=False)

    pd.set_option("display.width", 220)
    print("=== authoritative (20,000 resamples, seed 42) vs superseded (2,000, seed 0) ===")
    print(t[["threshold", "n", "model", "freq", "diff", "lo", "hi", "p", "ci_excludes_zero",
             "superseded_lo", "superseded_hi", "superseded_p", "superseded_ci_excludes_zero",
             "reproduces_frozen"]].to_string(index=False))
    print(f"\nall superseded values reproduce the frozen CSV: {bool(t.reproduces_frozen.all())}")
    flips = t[t.ci_excludes_zero != t.superseded_ci_excludes_zero]
    print(f"thresholds whose zero-crossing verdict changes: "
          f"{list(flips.threshold.values) if len(flips) else 'none'}")

    # ---- size effect with a proper paired interval, now that 1.01 arrays exist ----
    sm = pd.read_csv(SM, sep="\t")
    unr = fr[fr.threshold == 1.01].set_index("pdb")
    out = []
    for thr in (0.95, 0.85, 0.70):
        s = sm[sm.threshold == thr].copy()
        s["d"] = s.model_auprc - s.freq_auprc
        wide = s.pivot(index="pdb", columns="draw_seed", values="d")
        order = wide.index
        dbar = wide.mean(axis=1).values                       # size-matched principal
        dunr = (unr.loc[order, "model_auprc"] - unr.loc[order, "freq_auprc"]).values
        clu = s.drop_duplicates("pdb").set_index("pdb").loc[order, "clu70"].values
        d, lo, hi, p = boot(dunr, dbar, clu, **STD)
        out.append(dict(threshold=thr, unrestricted=round(float(np.nanmean(dunr)), 4),
                        size_matched=round(float(np.nanmean(dbar)), 4),
                        size_effect=round(d, 4), lo=round(lo, 4), hi=round(hi, 4), p=round(p, 4)))
    se = pd.DataFrame(out); se.to_csv(OUT/"size_effect_paired.tsv", sep="\t", index=False)
    print("\n=== size effect, PAIRED: unrestricted minus size-matched random, same complexes ===")
    print(se.to_string(index=False))
    print(f"\nwrote {OUT} and results/tables/similarity_curve_p84v2_authoritative.csv")


if __name__ == "__main__":
    main()
