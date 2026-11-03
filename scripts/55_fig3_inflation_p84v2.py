#!/usr/bin/env python3
"""Figure 3 (p84v2) - protocol-dependent performance inflation. p84v2 port of 25c_fig3_inflation_p82v1.py
(labels keyed by IMGT number + insertion code, tokens placed by position; results/review/CORRECTION_PLAN_P84V2.md).
Changes from 25c: inputs/outputs p84v2; panel B annotates the < 0.70 point (p = 0.052, CI touching zero)
instead of the p82v1 < 0.90 point; panel D y-range from the data; console cross-checks report p84v2 values.
(A) model vs frequency prior across Levenshtein-similarity thresholds
(B) model advantage with 95% cluster-bootstrap CIs (discrete point coloring only --
    no continuous significance shading, per JCIM_SCIENTIFIC_AND_FIGURE_AUDIT.md 9.2)
(C) docking-angle MAE, collapses under grouped evaluation
(D) docking-angle R-squared (newly persisted under p82v1) -- negative under both
    grouped protocols, positive under random
(E) training-size-matched control for the sweep (added 29 Oct 2026): at the three controlled
    thresholds, the similarity-filtered advantage, the five size-matched random draws and the
    principal size-matched summary with its cluster-bootstrap interval. Design frozen in
    results/review/SIZE_MATCHED_CONTROL_PLAN.md; estimands in SIZE_MATCHED_CONTROL_RESULT.md.
Sources: similarity_curve_p84v2.csv, angle_runs_p84v2.jsonl (both p84v2; the p82v1 source read
candidate set). No 28-position curve or point is plotted anywhere in this figure --
the 28-position pipeline is superseded, not a sensitivity comparator; see
results/review/P82V1_STAGE_A_NUMERICAL_AUDIT.md for the methodological audit of why.
"""
import json
import sys
sys.path.insert(0, ".")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from src.viz.style import apply_style, pstr, save_fig, print_validation_table, \
    MODEL_LABEL, SIM_LABEL, METHOD_COLORS, sig_marker, minus, REFERENCE_LINE_STYLE, \
    PANEL_LABEL_STYLE as PL, panel_letters

MODEL_C = METHOD_COLORS[MODEL_LABEL]
PRIOR_C = METHOD_COLORS["Pair-frequency prior"]
# C/D: predictions come from AngleNet, a separate regressor defined and trained in 16b (not
# TCR-IMAP), on 16b's ad hoc splits (not the frozen protocols) -> one neutral near-black,
# splits are named on the x-axis instead of colour-coded
ANGLE_C = METHOD_COLORS["AngleNet"]   # declared in style.py

apply_style()
PANELS = {}

# Authoritative uncertainty: every threshold recomputed at the project standard of 20,000
# resamples, seed 42, from the per-complex arrays recovered under a validation gate
# (results/review/SWEEP_UNCERTAINTY_RECONCILIATION.md). The frozen similarity_curve_p84v2.csv
# carries the run's own 2,000-resample seed-0 values and is retained as the superseded record;
# point estimates are identical in both.
SC_PATH = "results/tables/similarity_curve_p84v2_authoritative.csv"
SC_SUPERSEDED = "results/tables/similarity_curve_p84v2.csv"
ANG_PATH = "results/tables/angle_runs_p84v2.jsonl"
OUT_STEM = "figures/fig3_inflation_p84v2_v1"

# exact per-fold training sizes retained by the similarity filter, from
# results/review/control_estimands/similarity_diagnostics.tsv (regime = similarity-filtered)
TRSZ = {0.95: (91, 94, 96, 95, 94), 0.85: (86, 86, 88, 93, 86), 0.70: (72, 69, 76, 78, 73)}

sc = pd.read_csv(SC_PATH).sort_values("threshold")
_old = pd.read_csv(SC_SUPERSEDED).sort_values("threshold")
assert np.allclose(sc['diff'].values, _old['diff'].values), "point estimates must be unchanged"
assert np.allclose(sc.model.values, _old.model.values) and np.allclose(sc.freq.values, _old.freq.values)
ang = pd.DataFrame([json.loads(l) for l in open(ANG_PATH)])

assert (sc.n == 174).all(), "expected n=174 at every threshold"
assert (ang.n == 172).all(), "expected n=172 at every docking-angle row"

fig = plt.figure(figsize=(7.2, 7.9), constrained_layout=True)
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.05])

# --- A  model vs prior --------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
x = sc.threshold.values
ax.plot(x, sc.model, "o-", color=MODEL_C, lw=1.5, ms=4, label=MODEL_LABEL)
ax.plot(x, sc.freq, "s-", color=PRIOR_C, lw=1.3, ms=3.5, label="Pair-frequency prior")
ax.set_xlabel(f"Max train/test {SIM_LABEL} allowed")
ax.set_ylabel("Macro-AUPRC")
ax.set_xlim(1.06, 0.63)
ax.set_ylim(0.28, 0.68)
ax.set_xticks([1.01, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70])
ax.set_xticklabels(["none", "0.95", "0.90", "0.85", "0.80", "0.75", "0.70"], fontsize=6.2)
ax.legend(frameon=False, fontsize=6, loc="upper right", handlelength=1.2)
UNR = sc[sc.threshold > 1.0].iloc[0]
STR = sc[sc.threshold == 0.70].iloc[0]
ax.text(0.97, 0.60, minus("drop from none to 0.70: %s %+.4f, prior %+.4f"
        % (MODEL_LABEL, STR.model - UNR.model, STR.freq - UNR.freq)),
        transform=ax.transAxes, fontsize=6, ha="right", va="top",
        linespacing=1.4, color="0.25")
ax.text(0.03, 0.03, minus("n = 174 at every threshold"), transform=ax.transAxes,
        fontsize=6, color="0.35")
ax.set_title("Performance vs redundancy", fontsize=8, pad=3)
PANELS["A"] = [ax]

# --- B  advantage with CIs, discrete point coloring only, no shading ----
ax = fig.add_subplot(gs[0, 1])
# fill = the 95% CI excludes zero (paper-wide convention), NOT p < 0.05
ci_excl = ~((sc.lo <= 0) & (0 <= sc.hi))
ax.errorbar(x, sc["diff"], yerr=[sc["diff"] - sc.lo, sc.hi - sc["diff"]],
            fmt="none", ecolor="0.6", elinewidth=0.9, capsize=2.2)
for xi, di, ok in zip(x, sc["diff"], ci_excl):
    ax.plot(xi, di, "o", ms=5.2, zorder=3, **sig_marker(MODEL_C, bool(ok)))
ax.axhline(0, ls="--", lw=0.8, color="0.4")
ax.set_xlabel(f"Max train/test {SIM_LABEL} allowed")
ax.set_ylabel(f"{MODEL_LABEL} − prior (macro-AUPRC)")
ax.set_xlim(1.06, 0.63)
ax.set_xticks([1.01, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70])
ax.set_xticklabels(["none", "0.95", "0.90", "0.85", "0.80", "0.75", "0.70"], fontsize=6.2)
U = sc[sc.threshold > 1.0].iloc[0]
T = sc[sc.threshold == 0.70].iloc[0]
ax.text(1.005, U["diff"] - 0.030, minus("%+.4f\n%s" % (U["diff"], pstr(U.p))),
        fontsize=6, ha="left", va="top", color=MODEL_C, linespacing=1.35)
N90 = sc[sc.threshold == 0.90].iloc[0]
# p84v2: the borderline point is < 0.70 (CI touches zero); annotate it with its interval
ax.annotate(minus("%+.4f, %s\nCI [%+.4f, %+.4f]" % (T["diff"], pstr(T.p), T.lo, T.hi)),
            xy=(0.70, T["diff"]), xytext=(0.745, T["diff"] + 0.075),
            fontsize=5.8, ha="center", va="bottom", color="0.3", linespacing=1.3,
            arrowprops=dict(arrowstyle="-", lw=0.6, color="0.5", shrinkA=0, shrinkB=3))
ax.text(0.97, 0.93, minus("filled = 95% CI excludes zero"), transform=ax.transAxes,
        fontsize=6, color="0.35", ha="right", va="top")
ax.set_title("Advantage over a shortcut prior", fontsize=8, pad=3)
PANELS["B"] = [ax]

# --- C  docking-angle MAE ------------------------------------------------
ax = fig.add_subplot(gs[1, 0])
ORD = ["random", "allele", "trbv"]
LAB = {"random": "random", "allele": "allele\nholdout", "trbv": "TRBV\nholdout"}
base = ang.baseline_mae.iloc[0]
for i, s in enumerate(ORD):
    d = ang[ang.split_by == s]
    c = ANGLE_C
    ax.scatter(np.full(len(d), i) + np.linspace(-0.09, 0.09, len(d)),
               d.mae, s=24, color=c, zorder=3, edgecolor="none")
    ax.hlines(d.mae.mean(), i - 0.22, i + 0.22, color=c, lw=1.6, zorder=2)
    ax.text(i, 1.4, minus("r = %.2f to %.2f" %(d.pearson.min(), d.pearson.max())),
            ha="center", fontsize=6, color=c)
ax.axhline(base, ls="--", lw=1.0, color="0.35")
ax.text(-0.42, base + 0.35, minus("median baseline %.1f°" % base), fontsize=6,
        ha="left", va="bottom", color="0.35")
ax.set_xticks(range(3))
ax.set_xticklabels([LAB[s] for s in ORD], fontsize=6.5, linespacing=1.35)
ax.set_ylabel("Docking-angle MAE (°)")
ax.set_xlim(-0.5, 2.5)
ax.set_ylim(0, 18)
ax.text(0.03, 0.03, minus("n = 172 of 174 (3vxu, 3w0w excluded)"), transform=ax.transAxes,
        fontsize=6, color="0.35")
ax.set_title("Docking-angle regression", fontsize=8, pad=3)
PANELS["C"] = [ax]

# --- D  docking-angle R-squared (new: now persisted under p82v1) --------
ax = fig.add_subplot(gs[1, 1])
for i, s in enumerate(ORD):
    d = ang[ang.split_by == s]
    c = ANGLE_C
    ax.scatter(np.full(len(d), i) + np.linspace(-0.09, 0.09, len(d)),
               d.r2, s=24, color=c, zorder=3, edgecolor="none")
    ax.hlines(d.r2.mean(), i - 0.22, i + 0.22, color=c, lw=1.6, zorder=2)
ax.axhline(0, ls="--", lw=1.0, color="0.35")
ax.text(2.45, 0.06, minus("R² = 0\n(= mean predictor)"), fontsize=6, ha="right",
        va="bottom", color="0.35", linespacing=1.3)
ax.set_xticks(range(3))
ax.set_xticklabels([LAB[s] for s in ORD], fontsize=6.5, linespacing=1.35)
ax.set_ylabel("Docking-angle R²")
ax.set_xlim(-0.5, 2.5)
ax.set_ylim(min(-0.7, ang.r2.min() - 0.15), max(0.7, ang.r2.max() + 0.15))
ax.set_title("Docking-angle R²", fontsize=8, pad=3)
PANELS["D"] = [ax]

# --- E  training-size-matched control ------------------------------------
# Does the decline in B reflect homology control, or merely the smaller training sets that
# similarity filtering produces? Panel E answers that at the three controlled thresholds.
CE = pd.read_csv("results/review/control_estimands/control_estimands.tsv", sep="\t")
CE = CE.sort_values("threshold", ascending=False).reset_index(drop=True)
PD = pd.read_csv("results/review/size_matched_control/per_draw_metrics.tsv", sep="\t")
PD = PD[PD.draw_seed != "mean_over_draws"].copy()
PD["diff"] = PD["diff"].astype(float)

ax = fig.add_subplot(gs[2, :])
xs = np.arange(len(CE))
RAND_C = "0.30"

# unrestricted reference
ax.axhline(float(UNR["diff"]), zorder=1, **REFERENCE_LINE_STYLE)
ax.text(len(CE) - 0.5, float(UNR["diff"]) - 0.008,
        minus("unrestricted training (%d\u2013%d complexes per fold): %+.4f"
              % (139, 140, UNR["diff"])),
        fontsize=5.8, ha="right", va="top", color="0.35")

for i, r in CE.iterrows():
    thr = r.threshold
    # similarity-filtered estimate, with its interval under the common scheme
    ax.errorbar(i - 0.17, r.filtered_delta,
                yerr=[[r.filtered_delta - r.filtered_lo], [r.filtered_hi - r.filtered_delta]],
                fmt="none", ecolor="0.6", elinewidth=0.9, capsize=2.2, zorder=2)
    ax.plot(i - 0.17, r.filtered_delta, "s", ms=5.4, zorder=3,
            **sig_marker(MODEL_C, bool(r.filtered_lo > 0)))
    # five individual size-matched draws
    d = PD[np.isclose(PD.threshold, thr)]["diff"].values
    ax.scatter(np.full(len(d), i + 0.17) + np.linspace(-0.055, 0.055, len(d)), d,
               s=11, color=RAND_C, alpha=0.55, zorder=3, edgecolor="none")
    # principal size-matched summary
    ax.errorbar(i + 0.17, r.principal_mean_of_draw_differences,
                yerr=[[r.principal_mean_of_draw_differences - r.principal_lo],
                      [r.principal_hi - r.principal_mean_of_draw_differences]],
                fmt="none", ecolor="0.6", elinewidth=0.9, capsize=2.2, zorder=2)
    ax.plot(i + 0.17, r.principal_mean_of_draw_differences, "o", ms=5.6, zorder=4,
            **sig_marker(RAND_C, bool(r.principal_lo > 0)))
    # paired contrast between the two regimes
    ax.annotate("", xy=(i - 0.17, r.filtered_delta), xytext=(i + 0.17, r.principal_mean_of_draw_differences),
                arrowprops=dict(arrowstyle="-", lw=0.7, color="0.55", ls=(0, (2, 1.6))), zorder=1)
    ax.text(i, (r.filtered_delta + r.principal_mean_of_draw_differences) / 2,
            minus("paired\n%+.4f" % r.contrast_filtered_minus_sizematched),
            fontsize=5.6, ha="center", va="center", color="0.25", linespacing=1.25,
            bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none", alpha=0.82))

ax.axhline(0, ls="--", lw=0.8, color="0.4", zorder=1)
ax.set_xticks(xs)
ax.set_xticklabels([minus("%.2f\ntraining %d\u2013%d per fold"
                          % (r.threshold,
                             min(int(v) for v in TRSZ[r.threshold]),
                             max(int(v) for v in TRSZ[r.threshold])))
                    for _, r in CE.iterrows()], fontsize=6.2, linespacing=1.35)
ax.set_xlim(-0.55, len(CE) - 0.45)
ax.set_xlabel(f"Max train/test {SIM_LABEL} allowed (thresholds with a size-matched control)")
ax.set_ylabel(f"{MODEL_LABEL} \u2212 prior (macro-AUPRC)")
h = [plt.Line2D([], [], ls="none", marker="s", ms=5.4, mfc=MODEL_C, mec=MODEL_C,
                label="similarity-filtered"),
     plt.Line2D([], [], ls="none", marker="o", ms=5.6, mfc=RAND_C, mec=RAND_C,
                label="size-matched random (principal)"),
     plt.Line2D([], [], ls="none", marker="o", ms=3.4, mfc=RAND_C, mec="none", alpha=0.55,
                label="individual draws (5)")]
ax.legend(handles=h, frameon=False, fontsize=6, loc="upper left", handlelength=1.1,
          borderaxespad=0.3, ncol=3, columnspacing=1.4)
ax.set_title("Training-size-matched control: same training-set size, no similarity filter",
             fontsize=8, pad=3)
ax.set_ylim(-0.055, float(UNR["diff"]) + 0.105)
ax.text(0.5, 0.012, minus("error bars: 95% cluster bootstrap (20,000 resamples, seed 42); "
        "the five draws are shown as points, which are a spread and not an interval"),
        transform=ax.transAxes, fontsize=5.6, ha="center", va="bottom", color="0.4")
PANELS["E"] = [ax]

panel_letters(fig, PANELS)          # one fixed offset for every figure (style.py)
save_fig(fig, OUT_STEM)
plt.close(fig)

# ---------------- panel-to-source numerical validation table -------------
print("\n--- Panel A/B: similarity_curve_p84v2.csv (source of truth, read verbatim) ---")
print(sc[["threshold", "n", "model", "freq", "diff", "lo", "hi", "p"]].to_string(index=False))
print("\nDisplayed vs source, panel A/B extremes:")
print_validation_table(
    [dict(panel="A/B", threshold="unrestricted", n=int(UNR.n), model=UNR.model, freq=UNR.freq,
          diff=UNR["diff"], p_raw=UNR.p, p_displayed=pstr(UNR.p)),
     dict(panel="A/B", threshold="<0.90", n=int(sc[sc.threshold == 0.90].iloc[0].n),
          model=sc[sc.threshold == 0.90].iloc[0].model, freq=sc[sc.threshold == 0.90].iloc[0].freq,
          diff=sc[sc.threshold == 0.90].iloc[0]["diff"], p_raw=sc[sc.threshold == 0.90].iloc[0].p,
          p_displayed=pstr(sc[sc.threshold == 0.90].iloc[0].p)),
     dict(panel="A/B", threshold="<0.70", n=int(T.n), model=T.model, freq=T.freq,
          diff=T["diff"], p_raw=T.p, p_displayed=pstr(T.p))],
    [("panel", "panel", None), ("threshold", "threshold", None), ("n", "n", "%d"),
     ("model", "model_auprc", "%.4f"), ("freq", "prior_auprc", "%.4f"),
     ("diff", "diff", "%+.4f"), ("p_raw", "p_raw", "%.4f"), ("p_displayed", "p_displayed", None)])
print(f"\nNote (p84v2, authoritative 20,000 resamples seed 42): the < 0.70 point is {T['diff']:+.4f}, "
      f"p = {T.p:.3f}, 95% CI [{T.lo:+.4f}, {T.hi:+.4f}]. Under the superseded 2,000-resample seed-0 settings the same "
      "point estimate carried CI [-0.0001, +0.0543], p = 0.052, so the zero-crossing verdict at this threshold is not "
      "stable between two defensible bootstrap settings. Report the advantage as small and do NOT describe a "
      f"significance boundary. The < 0.90 point is p = {N90.p:.3f}, CI [{N90.lo:+.4f}, {N90.hi:+.4f}].")

print("\n--- Panel C/D: angle_runs_p84v2.jsonl (source of truth, read verbatim) ---")
print(ang[["split_by", "seed", "mae", "pearson", "r2", "n"]].to_string(index=False))
summary = ang.groupby("split_by").agg(
    mae_mean=("mae", "mean"), mae_min=("mae", "min"), mae_max=("mae", "max"),
    r_min=("pearson", "min"), r_max=("pearson", "max"),
    r2_min=("r2", "min"), r2_max=("r2", "max")).round(4)
print(summary.to_string())
print("baseline MAE %.3f (all rows identical: %s)" % (base, ang.baseline_mae.nunique() == 1))
print("\np84v2 docking-angle ranges (AngleNet, a separate regressor, not TCR-IMAP):")
for _s in ("random", "allele", "trbv"):
    _d = ang[ang.split_by == _s]
    print("  %-7s MAE [%.2f, %.2f]  r [%.3f, %.3f]  R2 [%+.3f, %+.3f]" % (_s, _d.mae.min(), _d.mae.max(), _d.pearson.min(), _d.pearson.max(), _d.r2.min(), _d.r2.max()))
