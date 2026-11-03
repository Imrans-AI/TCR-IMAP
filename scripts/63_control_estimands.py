#!/usr/bin/env python3
"""Size-matched control: the principal estimand, the ensemble estimand, and the direct paired
contrast against the similarity-filtered regime.

The committed plan (results/review/SIZE_MATCHED_CONTROL_PLAN.md) is AMBIGUOUS about the estimand:
section 3 says "Delta_rand -- this control, averaged over the five draws" (a mean of draw-specific
differences), while section 4 says the cluster bootstrap is "applied to the mean-over-draws
per-complex PREDICTIONS" (a prediction ensemble). Average precision is nonlinear, so these are not
the same quantity and they do not agree. The ambiguity is declared openly in the report; this
script computes both and labels the mean of draw-specific per-complex differences as PRINCIPAL,
with the prediction ensemble kept SECONDARY.

PRINCIPAL estimand, per complex c and threshold t:
    d_c(t) = (1/5) * sum over the 5 draws of [ AP_model(c, draw) - AP_prior(c, draw) ]
bootstrapped over clu70 clusters (20,000 resamples, seed 42). It is CONDITIONAL on the five
realized training draws; the draw-to-draw range is reported separately and never merged into it.

DIRECT PAIRED CONTRAST, per complex c and threshold t:
    g_c(t) = [ AP_model_filtered(c) - AP_prior_filtered(c) ] - d_c(t)
on the same 174 test complexes, bootstrapped with the same clusters and resamples, so the two
regimes are compared under one resampling scheme rather than by differencing two separate
intervals.

usage: python scripts/63_control_estimands.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, ".")
ROOT = Path(__file__).resolve().parent.parent
SM = ROOT/"results/review/size_matched_control"
FR = ROOT/"results/review/filtered_rerun"
OUT = ROOT/"results/review/control_estimands"; OUT.mkdir(parents=True, exist_ok=True)
THRESHOLDS = [0.95, 0.85, 0.70]
N_BOOT, SEED = 20000, 42


def boot(vals, clusters, n_boot=N_BOOT, seed=SEED):
    """Cluster bootstrap of the mean of a per-complex quantity. Returns
    (mean, lo, hi, p_exact, p_hit_floor). p is two-sided about zero; the exact value is kept at
    full precision and the floor is reported separately from it."""
    rng = np.random.default_rng(seed)
    ok = ~np.isnan(vals)
    v, c = np.asarray(vals)[ok], np.asarray(clusters)[ok]
    uc = np.unique(c)
    idx_by = {u: np.where(c == u)[0] for u in uc}
    obs = v.mean()
    d = np.empty(n_boot)
    for i in range(n_boot):
        pick = rng.choice(uc, len(uc), replace=True)
        d[i] = v[np.concatenate([idx_by[u] for u in pick])].mean()
    lo, hi = np.percentile(d, [2.5, 97.5])
    n_le, n_ge = int((d <= 0).sum()), int((d >= 0).sum())
    p_raw = 2*min(n_le, n_ge)/n_boot
    return dict(mean=obs, lo=lo, hi=hi, n_le=n_le, n_ge=n_ge,
                p_raw=p_raw, p_floor=1.0/n_boot, p_reported=max(p_raw, 1.0/n_boot),
                hit_floor=bool(p_raw < 1.0/n_boot))


def main():
    sm = pd.read_csv(SM/"per_complex_metrics.tsv", sep="\t")
    fr = pd.read_csv(FR/"filtered_per_complex.tsv", sep="\t")
    draws = sorted(sm.draw_seed.unique())
    assert len(draws) == 5, draws
    per_draw = pd.read_csv(SM/"per_draw_metrics.tsv", sep="\t")

    rows, pc_rows = [], []
    for thr in THRESHOLDS:
        a = sm[sm.threshold == thr].copy()
        b = fr[fr.threshold == thr].set_index("pdb")
        a["d"] = a.model_auprc - a.freq_auprc
        wide = a.pivot(index="pdb", columns="draw_seed", values="d")
        assert wide.shape[1] == 5
        order = wide.index
        clusters = a.drop_duplicates("pdb").set_index("pdb").loc[order, "clu70"].values
        d_bar = wide.mean(axis=1).values                      # PRINCIPAL per-complex quantity
        d_filt = (b.loc[order, "model_auprc"] - b.loc[order, "freq_auprc"]).values
        g = d_filt - d_bar                                    # DIRECT paired contrast

        P = boot(d_bar, clusters)                             # principal
        G = boot(g, clusters)                                 # direct contrast
        F = boot(d_filt, clusters)                            # filtered regime, same scheme

        gr = per_draw[(per_draw.threshold == thr) & (per_draw.draw_seed != "mean_over_draws")]
        ens = per_draw[(per_draw.threshold == thr) & (per_draw.draw_seed == "mean_over_draws")].iloc[0]

        rows.append(dict(
            threshold=thr,
            # principal
            principal_mean_of_draw_differences=round(P["mean"], 4),
            principal_lo=round(P["lo"], 4), principal_hi=round(P["hi"], 4),
            principal_p_raw=P["p_raw"], principal_p_reported=P["p_reported"],
            principal_hit_floor=P["hit_floor"],
            # draw spread, kept separate
            draw_min=round(gr["diff"].min(), 4), draw_max=round(gr["diff"].max(), 4),
            draw_range=round(gr["diff"].max()-gr["diff"].min(), 4),
            # filtered regime under the same scheme
            filtered_delta=round(F["mean"], 4), filtered_lo=round(F["lo"], 4),
            filtered_hi=round(F["hi"], 4),
            # the direct paired contrast
            contrast_filtered_minus_sizematched=round(G["mean"], 4),
            contrast_lo=round(G["lo"], 4), contrast_hi=round(G["hi"], 4),
            contrast_p_raw=G["p_raw"], contrast_p_reported=G["p_reported"],
            contrast_hit_floor=G["hit_floor"],
            # secondary: prediction ensemble
            secondary_ensemble_delta=ens["diff"], secondary_ensemble_lo=ens.lo,
            secondary_ensemble_hi=ens.hi))
        for k, p in enumerate(order):
            pc_rows.append(dict(threshold=thr, pdb=p, clu70=clusters[k],
                                mean_draw_difference=d_bar[k], filtered_difference=d_filt[k],
                                contrast=g[k]))

    t = pd.DataFrame(rows)
    t.to_csv(OUT/"control_estimands.tsv", sep="\t", index=False)
    pd.DataFrame(pc_rows).to_csv(OUT/"per_complex_estimands.tsv", sep="\t", index=False)

    pd.set_option("display.width", 200)
    print("=== PRINCIPAL: mean over draws of per-complex (AP_model - AP_prior), bootstrapped by clu70 ===")
    print(t[["threshold", "principal_mean_of_draw_differences", "principal_lo", "principal_hi",
             "principal_p_raw", "principal_p_reported", "principal_hit_floor",
             "draw_min", "draw_max", "draw_range"]].to_string(index=False))
    print("\n=== SECONDARY: prediction-ensemble (not the principal estimand) ===")
    print(t[["threshold", "secondary_ensemble_delta", "secondary_ensemble_lo",
             "secondary_ensemble_hi"]].to_string(index=False))
    print("\n=== DIRECT PAIRED CONTRAST: filtered minus size-matched, same complexes and resamples ===")
    print(t[["threshold", "filtered_delta", "filtered_lo", "filtered_hi",
             "contrast_filtered_minus_sizematched", "contrast_lo", "contrast_hi",
             "contrast_p_raw", "contrast_p_reported", "contrast_hit_floor"]].to_string(index=False))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
