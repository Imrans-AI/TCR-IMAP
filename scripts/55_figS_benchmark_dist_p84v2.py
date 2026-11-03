#!/usr/bin/env python3
"""SI figure (p84v2) - per-complex per-residue AUPRC distributions behind Figure 4A.

Built from 37e_figS_benchmark_dist_p82v1.py; inputs are the p84v2 per-complex arrays.

One panel per protocol (random, double, allele n = 64, temporal, tcr70), five methods per panel:
TCR-IMAP, Struct-TCS, Seq-TCS, position prior, NN transfer. Box = median and interquartile range,
whiskers 1.5 x IQR; every complex drawn as a point (horizontal jitter for display only, fixed
seed). Box outlines and points in METHOD_COLORS.

Sources: results/tables/percomplex_{random,double,allele,temporal,tcr70}_p84v2.npz (written by
scripts/53_benchmark_p84v2.py), checked against results/tables/benchmark_p84v2.csv: every
cell's mean over non-missing complexes must equal the recorded macro-AUPRC to 4 decimals.
Nothing is recomputed: these are the stored per-complex values. A method-protocol cell whose
array is missing would be omitted and reported (none is missing). Complexes with an undefined
AUPRC (no positive or no negative residue) are stored as NaN and dropped, and their count is
reported per cell.
No statistic is tested in this figure, so no marker carries the fill = significance meaning.
"""
import sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, "src")
from viz.style import (apply_style, save_fig, print_validation_table, PAGE_WIDTH_IN,
                       MODEL_LABEL, METHOD_COLORS, minus)
from viz.figstats import pdf_font_check

apply_style()
T = "results/tables/"
OUT_STEM = "figures/figS_benchmark_distributions_p84v2_v1"
PROTOS = ["random", "double", "allele", "temporal", "tcr70"]       # light -> dark strictness
KEY2LAB = {"Ours": MODEL_LABEL, "Struct-TCS": "Struct-TCS", "Seq-TCS": "Seq-TCS",
           "PosFreq": "Position prior", "NN": "NN transfer"}
METHODS = list(KEY2LAB)
C = {k: METHOD_COLORS[v] for k, v in KEY2LAB.items()}

bench = pd.read_csv(T + "benchmark_p84v2.csv").set_index("protocol")
data, omitted, rows = {}, [], []
for p in PROTOS:
    z = np.load(T + f"percomplex_{p}_p84v2.npz", allow_pickle=True)
    assert len(z["pdb"]) == int(bench.loc[p, "n"]), p
    for k in METHODS:
        key = f"{k}_auprc"
        if key not in z.files:
            omitted.append((p, k)); continue
        v = np.asarray(z[key], float)
        ok = v[~np.isnan(v)]
        data[(p, k)] = ok
        mean = float(ok.mean())
        assert round(mean, 4) == bench.loc[p, k], (p, k, mean, bench.loc[p, k])
        q1, med, q3 = np.percentile(ok, [25, 50, 75])
        rows.append(dict(protocol=p, method=KEY2LAB[k], n=len(v), nan=int(np.isnan(v).sum()),
                         mean=mean, rec=bench.loc[p, k], median=med, q1=q1, q3=q3,
                         lo=float(ok.min()), hi=float(ok.max())))
assert len(pcx_n := {p: int(bench.loc[p, "n"]) for p in PROTOS}) == 5 and pcx_n["allele"] == 64

FW, FH = PAGE_WIDTH_IN, 2.75
fig = plt.figure(figsize=(FW, FH))
L0, W, GAP, B, H = 0.55, 1.22, 0.10, 0.62, 1.70
rng = np.random.default_rng(0)        # jitter only
axes = []
for i, p in enumerate(PROTOS):
    ax = fig.add_axes([(L0 + i * (W + GAP)) / FW, B / FH, W / FW, H / FH], sharey=axes[0] if axes else None)
    axes.append(ax)
    for j, k in enumerate(METHODS):
        if (p, k) not in data:
            continue
        v = data[(p, k)]
        ax.scatter(j + rng.uniform(-0.22, 0.22, len(v)), v, s=2.2, color=C[k], alpha=0.45,
                   edgecolor="none", zorder=2, rasterized=False)
        ax.boxplot([v], positions=[j], widths=0.62, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor="none", edgecolor=C[k], lw=0.9),
                   medianprops=dict(color=C[k], lw=1.4), whiskerprops=dict(color=C[k], lw=0.8),
                   capprops=dict(color=C[k], lw=0.8), zorder=3)
    ax.set_xticks([]); ax.set_xlim(-0.6, len(METHODS) - 0.4)
    ax.set_title(minus(f"{p}\nn = {pcx_n[p]}" + (" (primary)" if p == "tcr70" else "")),
                 fontsize=7, pad=3, linespacing=1.2)
    if i == 0:
        ax.set_ylabel("Per-complex per-residue AUPRC")
    else:
        plt.setp(ax.get_yticklabels(), visible=False)
axes[0].set_ylim(0, 1.02)
fig.legend(handles=[Line2D([], [], marker="s", ls="", ms=5.5, mfc=C[k], mec="none",
                           label=KEY2LAB[k]) for k in METHODS],
           loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=5, frameon=False, fontsize=6.3,
           handletextpad=0.3, columnspacing=1.3)
fig.text(0.99, 0.18 / FH, minus("box: median and IQR, whiskers 1.5 × IQR; points: complexes"),
         ha="right", va="bottom", fontsize=5.6, color="0.40")

png, pdf = save_fig(fig, OUT_STEM)
plt.close(fig)

print("\n--- per-complex per-residue AUPRC (percomplex_*_p84v2.npz) vs benchmark_p84v2.csv ---")
print_validation_table(rows, [("protocol", "protocol", None), ("method", "method", None),
                              ("n", "n", "%d"), ("nan", "NaN dropped", "%d"),
                              ("mean", "mean", "%.4f"), ("rec", "recorded macro", "%.4f"),
                              ("median", "median", "%.4f"), ("q1", "Q1", "%.4f"), ("q3", "Q3", "%.4f"),
                              ("lo", "min", "%.4f"), ("hi", "max", "%.4f")])
print(f"\nomitted method-protocol cells (no per-complex array): {omitted if omitted else 'none'}")
# tcr70 means under the cluster-disjoint validation split (CORRECTION_PLAN_P84V2.md, "Per-residue validation
# split"). Under the superseded random split they were 0.7446 / 0.6942 / 0.6548 / 0.7196; the position prior is
# unchanged by construction, since it is fitted on training and validation together.
EXP = {"Ours": 0.7281, "Struct-TCS": 0.6969, "Seq-TCS": 0.6481, "PosFreq": 0.7196}
print("benchmark_p84v2.csv cross-check (per-residue, tcr70 means):")
for k, e in EXP.items():
    m = float(data[("tcr70", k)].mean())
    print(f"  {'OK      ' if round(m, 4) == e else 'MISMATCH'} {KEY2LAB[k]} {m:.4f} (recorded {e})")
    assert round(m, 4) == e
n3, arial, bases = pdf_font_check(pdf)
print(f"\nPDF font check {pdf}: Type3 = {n3}, ArialMT present = {arial}, BaseFonts = {bases}")
assert n3 == 0 and arial
