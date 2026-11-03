#!/usr/bin/env python3
"""SI figure S9 (p84v2) - the identity-threshold sweep of Figure 3A/B with its composition.

A  TCR-IMAP - pair-frequency prior (pairwise macro-AUPRC) per threshold, 95 % CI and p: Figure 3B's values
B  complexes per fold: test (unchanged at every threshold) and the training complexes each fold kept
C  train/test similarity achieved: per test complex, the maximum similarity to the training complexes its fold kept;
   the mean and the maximum over the 174 test complexes; with the threshold itself and, for comparison, the existing
   mean_max_sim column of similarity_curve_p84v2.csv, which is an upper bound (see below), not the achieved value

Sources (nothing is trained, predicted or re-scored):
- results/tables/similarity_curve_p84v2.csv: Delta, CI, p (recorded, 2,000 resamples, seed 0; the sweep stores no
  per-complex arrays) and n
- results/tables/similarity_curve_composition_p84v2.csv and _folds_p84v2.csv, written by
  scripts/57_sweep_composition_p84v2.py from the committed similarity matrix and the sweep's frozen fold seed
  (authorised 26 October 2026)
The mean_max_sim column is, for each test complex, the maximum similarity to ANY other complex of the 174 below the
threshold - including complexes in the same test fold and training complexes the filter removed - averaged over test
complexes. It is therefore an upper bound on the achieved mean.
Fill = 95 % CI excludes zero (panel A only). B and C are descriptive.
"""
import sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, "src")
from viz.style import (apply_style, save_fig, print_validation_table, pstr, sig_marker, signed, minus, panel_letters,
                       PAGE_WIDTH_IN, MODEL_LABEL, METHOD_COLORS, SIM_LABEL, REFERENCE_LINE_STYLE)
from viz.figstats import ci_excludes_zero, pdf_font_check

apply_style()
T = "results/tables/"
OUT_STEM = "figures/figS_sweep_composition_p84v2_v1"
TXT = "0.30"
RD = dict(float_precision="round_trip")
XT = [1.01, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70]
XL = ["none", "0.95", "0.90", "0.85", "0.80", "0.75", "0.70"]
VAL = []
def val(panel, what, shown, src, cell):
    VAL.append(dict(panel=panel, what=what, shown=shown, src=src, cell=cell))

sc = pd.read_csv(T + "similarity_curve_p84v2.csv", **RD).sort_values("threshold", ascending=False).reset_index(drop=True)
co = pd.read_csv(T + "similarity_curve_composition_p84v2.csv", **RD)
fo = pd.read_csv(T + "similarity_curve_composition_folds_p84v2.csv", **RD)
assert np.allclose(sc.threshold, XT) and np.allclose(co.threshold, XT)
assert (sc.n == 174).all() and (co.n_test == 174).all() and (co.folds_used == 5).all()
assert np.allclose(co.mean_max_sim_column.round(4), sc.mean_max_sim)          # the column, reproduced by 57_sweep
assert (co.achieved_max_sim_mean < co.mean_max_sim_column).all()             # upper bound, strictly above here
assert (co.achieved_max_sim_max.values[1:] < np.array(XT[1:])).all()          # the cap holds
assert (fo.groupby("threshold").n_test.apply(sorted).apply(lambda v: v == [34, 35, 35, 35, 35])).all()
print("inputs consistent: n test 174 at all 7 thresholds; achieved mean below the mean_max_sim column everywhere; "
      "achieved maximum below every finite threshold")

FW, FH = PAGE_WIDTH_IN, 3.05
fig = plt.figure(figsize=(FW, FH))
def axin(l, b, w, h): return fig.add_axes([l / FW, b / FH, w / FW, h / FH])
def xaxis(ax):
    ax.set_xlim(1.06, 0.63); ax.set_xticks(XT); ax.set_xticklabels(XL, fontsize=6.2)
    ax.set_xlabel(f"Max train/test {SIM_LABEL}\nallowed", fontsize=7, linespacing=1.2)
CMOD = METHOD_COLORS[MODEL_LABEL]

# --- A
ax = axA = axin(0.72, 0.78, 1.62, 1.85)
ax.axhline(0, ls="--", lw=0.8, color="0.4")
for _, r in sc.iterrows():
    sig = ci_excludes_zero(r.lo, r.hi)
    ax.plot([r.threshold] * 2, [r.lo, r.hi], color="0.6", lw=0.9)
    ax.plot(r.threshold, r["diff"], "o", ms=4.0, zorder=3, **sig_marker(CMOD, sig))
    val("A", f"threshold {r.threshold:g}", f"{signed(r['diff'], 4)} [{signed(r.lo, 4)}, {signed(r.hi, 4)}] {pstr(r.p)} "
        f"({'filled' if sig else 'hollow'})", "similarity_curve_p84v2.csv", f"threshold={r.threshold:g}: diff, lo, hi, p")
r70 = sc.iloc[-1]
ax.annotate(minus(f"{signed(r70['diff'], 4)}, {pstr(r70.p)}"), xy=(0.70, r70["diff"]), xytext=(0.80, 0.20),
            fontsize=5.8, color=TXT, ha="center", arrowprops=dict(arrowstyle="-", color="0.5", lw=0.6))
xaxis(ax)
ax.set_ylabel(f"{MODEL_LABEL} − prior\n(pairwise macro-AUPRC)", linespacing=1.2)
ax.set_title("Advantage over the prior", fontsize=8, pad=4)

# --- B
ax = axB = axin(3.05, 0.78, 1.65, 1.85)
for thr, g in fo.groupby("threshold"):
    ax.plot([thr] * len(g), g.n_train_kept, "x", ms=3.6, mew=0.8, color="k", zorder=3)
    ax.plot([thr] * len(g), g.n_test, "_", ms=6, mew=1.0, color="0.15", zorder=3)
ax.plot(co.threshold, co.n_train_mean, "-", color="k", lw=1.0)
ax.plot(sc.threshold, [34.8] * len(sc), ls=(0, (3, 2)), color="0.15", lw=0.9)
for _, r in co.iterrows():
    val("B", f"threshold {r.threshold:g} training kept per fold", f"{int(r.n_train_min)}–{int(r.n_train_max)} (mean {r.n_train_mean:.1f})",
        "similarity_curve_composition_p84v2.csv", f"threshold={r.threshold:g}: n_train_min/max/mean")
ax.text(0.97, 0.10, "test: 34–35 per fold, 174 in all,\nat every threshold", transform=ax.transAxes, ha="right",
        va="bottom", fontsize=5.8, color=TXT, linespacing=1.2)
r1, r7 = co.iloc[0], co.iloc[-1]
ax.text(0.97, 0.97, minus(f"training kept per fold:\n{int(r1.n_train_min)}–{int(r1.n_train_max)} (none) → "
                          f"{int(r7.n_train_min)}–{int(r7.n_train_max)} (0.70)"),
        transform=ax.transAxes, ha="right", va="top", fontsize=5.8, color=TXT, linespacing=1.2)
val("B", "test per fold", "34–35, 174 in all", "similarity_curve_composition_folds_p84v2.csv; similarity_curve_p84v2.csv n", "n_test")
xaxis(ax); ax.set_ylim(0, 165)
ax.set_ylabel("Complexes per fold")
ax.set_title("Fold composition", fontsize=8, pad=4)
ax.legend(handles=[Line2D([], [], marker="x", ls="-", color="k", ms=3.6, mew=0.8, lw=1.0, label="training kept (folds, mean)"),
                   Line2D([], [], marker="_", ls=(0, (3, 2)), color="0.15", ms=6, lw=0.9, label="test")],
          loc="lower left", bbox_to_anchor=(0.0, 0.26), frameon=False, fontsize=5.8, handlelength=2.0)

# --- C
ax = axC = axin(5.45, 0.78, 1.65, 1.85)
fin = co.threshold < 1.0
ax.plot(co.threshold[fin], co.threshold[fin], **REFERENCE_LINE_STYLE)
ax.plot(co.threshold, co.achieved_max_sim_max, "^-", color="k", ms=3.4, lw=1.0)
ax.plot(co.threshold, co.achieved_max_sim_mean, "o-", color="k", ms=3.0, lw=1.0)
ax.plot(co.threshold, co.mean_max_sim_column, ls=(0, (4, 1.5, 1, 1.5)), color="0.45", lw=1.0)
for _, r in co.iterrows():
    val("C", f"threshold {r.threshold:g}", f"achieved max {r.achieved_max_sim_max:.4f}, mean {r.achieved_max_sim_mean:.4f}; "
        f"column {r.mean_max_sim_column:.4f}", "similarity_curve_composition_p84v2.csv",
        f"threshold={r.threshold:g}: achieved_max_sim_max, achieved_max_sim_mean, mean_max_sim_column")
xaxis(ax); ax.set_ylim(0.45, 1.02)
ax.set_ylabel("Train/test similarity", fontsize=7)
ax.set_title("Similarity achieved", fontsize=8, pad=4)
ax.legend(handles=[Line2D([], [], marker="^", color="k", ms=3.4, lw=1.0, label="achieved, maximum"),
                   Line2D([], [], marker="o", color="k", ms=3.0, lw=1.0, label="achieved, mean"),
                   Line2D([], [], **REFERENCE_LINE_STYLE, label="threshold"),
                   Line2D([], [], ls=(0, (4, 1.5, 1, 1.5)), color="0.45", lw=1.0, label="mean_max_sim (upper bound)")],
          loc="lower left", frameon=False, fontsize=5.6, handlelength=2.2)

fig.text(0.99, 0.008, "A: recorded values, clu70 cluster bootstrap, 2,000 resamples (Figure 3B). B, C: descriptive",
         ha="right", va="bottom", fontsize=5.4, color="0.45")
panel_letters(fig, {"A": [axA], "B": [axB], "C": [axC]})
png, pdf = save_fig(fig, OUT_STEM)
plt.close(fig)

print("\n--- every drawn number and its source ---")
print_validation_table(VAL, [("panel", "panel", None), ("what", "value", None), ("shown", "drawn", None),
                             ("src", "source", None), ("cell", "cell", None)])
n3, arial, bases = pdf_font_check(pdf)
print(f"\nPDF font check {pdf}: Type3 = {n3}, ArialMT present = {arial}, BaseFonts = {bases}")
assert n3 == 0 and arial
