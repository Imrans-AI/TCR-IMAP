#!/usr/bin/env python3
"""SI figure S5 (p84v2) - cluster-macro sensitivity of the three primary tcr70 comparisons.

Each comparison under complex-macro weighting (every test complex weighted equally; the predeclared primary) and under
cluster-macro weighting (each of the 69 clu70 clusters weighted equally, complexes averaged within their cluster first):
per-residue TCR-IMAP vs Struct-TCS, per-residue TCR-IMAP vs the position prior, pairwise TCR-IMAP vs the pair-frequency
prior. Difference with 95 % CI and p; the two methods' macro-AUPRC under each weighting.

Sources (nothing is trained, predicted or re-scored):
- results/tables/percomplex_tcr70_p84v2.npz (per-residue, cluster-disjoint validation split)
- results/tables/percomplex_tcr70_paired_p84v2.tsv (pairwise; checked against pairwise_ours_p84v2.npz)
- results/review/figures_p84v2/stats_20k_seed42.tsv: every drawn difference, CI and p is recomputed here with
  figstats.pooled_cluster_bootstrap (complex-macro) and figstats.cluster_macro_bootstrap (cluster-macro) at 20,000
  resamples, seed 42, and asserted equal to that table
- results/tables/cluster_macro_p84v2.csv: complex- and cluster-macro means (4 dp) and the recorded 2,000-resample
  cluster-macro rows, both reproduced here from the same arrays
Fill = 95 % CI excludes zero. Marker = weighting: circle complex-macro, square cluster-macro (as Figure 4E).
"""
import sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, "src")
from viz.style import (apply_style, save_fig, print_validation_table, pstr, sig_marker, signed, minus, panel_letters,
                       PAGE_WIDTH_IN, MODEL_LABEL, METHOD_COLORS)
from viz.figstats import (pooled_cluster_bootstrap, cluster_macro_bootstrap, ci_excludes_zero, pdf_font_check,
                          N_BOOT, SEED)

apply_style()
T = "results/tables/"
OUT_STEM = "figures/figS_cluster_macro_p84v2_v1"
TXT = "0.30"
RD = dict(float_precision="round_trip")
VAL = []
def val(what, shown, src, cell):
    VAL.append(dict(what=what, shown=shown, src=src, cell=cell))

# ------------------------------------------------------------------ inputs
z = np.load(T + "percomplex_tcr70_p84v2.npz", allow_pickle=True)
pr = pd.read_csv(T + "percomplex_tcr70_paired_p84v2.tsv", sep="\t", **RD)
po = np.load(T + "pairwise_ours_p84v2.npz", allow_pickle=True)
assert (po["pdb"] == pr.pdb.values).all() and np.allclose(po["auprc"], pr.grip.values) and (po["clu70"] == pr.clu70.values).all()
md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
assert (md.loc[z["pdb"], "clu70"].values == z["clu70"]).all() and md.clu70.nunique() == 69

ROWS = [  # key, label, a, b, clusters, stats_20k comparison names (complex, cluster), cluster_macro_p84v2 rows, colour
    ("struct", "per-residue\nvs Struct-TCS", z["Ours_auprc"], z["Struct-TCS_auprc"], z["clu70"],
     ("tcr70: TCR-IMAP vs Struct-TCS", "per-residue tcr70: TCR-IMAP vs Struct-TCS"),
     ("per-residue", "Ours", "Struct-TCS", "Ours_minus_Struct-TCS"), METHOD_COLORS["Struct-TCS"]),
    ("posfreq", "per-residue\nvs position prior", z["Ours_auprc"], z["PosFreq_auprc"], z["clu70"],
     ("tcr70: TCR-IMAP vs PosFreq", "per-residue tcr70: TCR-IMAP vs PosFreq"),
     ("per-residue", "Ours", "PosFreq", "Ours_minus_PosFreq"), METHOD_COLORS["Position prior"]),
    ("pair", "pairwise\nvs pair-frequency prior", pr.grip.values, pr.prior.values, pr.clu70.values,
     ("tcr70: TCR-IMAP vs pair-frequency prior", "pairwise tcr70: TCR-IMAP vs pair-frequency prior"),
     ("pairwise", "model(GRIP)", "freq(prior)", "model_minus_freq"), METHOD_COLORS["Pair-frequency prior"]),
]
S20 = pd.read_csv("results/review/figures_p84v2/stats_20k_seed42.tsv", sep="\t", **RD)
CM = pd.read_csv(T + "cluster_macro_p84v2.csv", **RD)
res = {}
for key, lab, a, b, cl, (n_cx, n_cl), (task, ma, mb, mdiff), col in ROWS:
    cx = pooled_cluster_bootstrap(a, b, cl, N_BOOT, SEED)
    cm = cluster_macro_bootstrap(a, b, cl, N_BOOT, SEED)
    for r, grp, name in ((cx, "per-residue" if task == "per-residue" else "pairwise", n_cx), (cm, "cluster-macro", n_cl)):
        rec = S20[(S20.group == grp) & (S20.comparison == name)].iloc[0]
        for f in ("diff", "lo", "hi", "p"):
            assert abs(r[f] - rec[f]) < 1e-12, (key, grp, f, r[f], rec[f])
    # means and the recorded 2k cluster-macro row
    c = CM[(CM.task == task) & (CM.protocol == "tcr70")].set_index("model")
    assert round(float(np.nanmean(a)), 4) == c.loc[ma, "complex_macro"] and round(float(np.nanmean(b)), 4) == c.loc[mb, "complex_macro"]
    assert round(cm["macro_a"], 4) == c.loc[ma, "cluster_macro"] and round(cm["macro_b"], 4) == c.loc[mb, "cluster_macro"], key
    r2 = cluster_macro_bootstrap(a, b, cl, 2000, 0)
    rr = c.loc[mdiff]
    assert (round(r2["diff"], 4), round(r2["lo"], 4), round(r2["hi"], 4), round(r2["p"], 4)) == \
        (rr.paired_diff_vs_GRIP_cluster_macro, rr.ci_lo_cluster_macro, rr.ci_hi_cluster_macro, rr.p_cluster_macro), key
    res[key] = dict(cx=cx, cm=cm, means=dict(cx=(float(np.nanmean(a)), float(np.nanmean(b))), cm=(cm["macro_a"], cm["macro_b"])))
print(f"3 comparisons x 2 weightings recomputed at {N_BOOT:,} resamples, seed {SEED}: all equal stats_20k_seed42.tsv; "
      "means and the recorded 2,000-resample cluster-macro rows reproduce cluster_macro_p84v2.csv")
agree = all(ci_excludes_zero(res[k]["cx"]["lo"], res[k]["cx"]["hi"]) == ci_excludes_zero(res[k]["cm"]["lo"], res[k]["cm"]["hi"])
            for k in res)
assert agree, "complex- and cluster-macro verdicts disagree - the caption's 'now agree' would be false"
print("verdict (CI excludes zero) agrees between the two weightings in all three comparisons")

# ------------------------------------------------------------------ figure
FW, FH = PAGE_WIDTH_IN, 2.75
fig = plt.figure(figsize=(FW, FH))
L, W = 1.55, 1.85
ax = fig.add_axes([L / FW, 0.52 / FH, W / FW, 1.80 / FH])
XD, XP, XM = 1 + 0.12 / W, 1 + 1.45 / W, 1 + 2.12 / W
WT = {"cx": ("complex-macro (primary)", "o"), "cm": ("cluster-macro", "s")}
yy = np.arange(len(ROWS))[::-1].astype(float)
for (key, lab, *_rest, col), y in zip(ROWS, yy):
    for dy, w in ((+0.17, "cx"), (-0.17, "cm")):
        r = res[key][w]; ma, mb = res[key]["means"][w]
        ax.plot([r["lo"], r["hi"]], [y + dy, y + dy], color=col, lw=1.0, solid_capstyle="butt")
        ax.plot(r["diff"], y + dy, WT[w][1], ms=4.4, zorder=3, **sig_marker(col, ci_excludes_zero(r["lo"], r["hi"])))
        ax.text(XD, y + dy, minus(f"{signed(r['diff'], 4)} [{signed(r['lo'], 4)}, {signed(r['hi'], 4)}]"),
                transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.8, color=TXT)
        ax.text(XP, y + dy, pstr(r["p"]), transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.8, color=TXT)
        ax.text(XM, y + dy, f"{ma:.4f} vs {mb:.4f}", transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=5.8, color=TXT)
        val(f"{lab.replace(chr(10), ' ')}, {WT[w][0]}",
            f"{signed(r['diff'], 4)} [{signed(r['lo'], 4)}, {signed(r['hi'], 4)}] {pstr(r['p'])}; {ma:.4f} vs {mb:.4f} "
            f"({'filled' if ci_excludes_zero(r['lo'], r['hi']) else 'hollow'})", "stats_20k_seed42.tsv; cluster_macro_p84v2.csv",
            f"{'cluster-macro' if w == 'cm' else 'per-residue/pairwise'} row; {'cluster_macro' if w == 'cm' else 'complex_macro'} means")
ax.axvline(0, ls="--", lw=0.8, color="0.4")
ax.set_yticks(yy); ax.set_yticklabels([r[1] for r in ROWS], fontsize=6.5, linespacing=1.15)
ax.set_ylim(-0.6, len(ROWS) - 0.4)
ax.set_xlabel(f"{MODEL_LABEL} − comparator (macro-AUPRC, tcr70)", fontsize=7)
for x, h in ((XD, "Δ [95% CI]"), (XP, "p"), (XM, f"{MODEL_LABEL} vs comparator")):
    ax.text(x, 1.0, h, transform=ax.transAxes, ha="left", va="bottom", fontsize=6.2, color="0.2")
ax.legend(handles=[Line2D([], [], marker=m, ls="", mfc="0.3", mec="0.3", ms=4.2, label=l) for l, m in WT.values()],
          loc="upper left", bbox_to_anchor=(-0.78, -0.17), ncol=2, frameon=False, fontsize=6.3, handletextpad=0.3,
          columnspacing=1.2)
ax.set_title("Cluster weighting: tcr70 primary comparisons", fontsize=8, pad=12)
fig.text(0.99, 0.01, minus(f"clu70 cluster bootstrap, 69 clusters, {N_BOOT:,} resamples, seed {SEED}, recomputed from stored "
                           "per-complex arrays"), ha="right", va="bottom", fontsize=5.4, color="0.45")
png, pdf = save_fig(fig, OUT_STEM)
plt.close(fig)

print("\n--- every drawn number and its source ---")
print_validation_table(VAL, [("what", "row", None), ("shown", "drawn", None), ("src", "source", None), ("cell", "cell", None)])
n3, arial, bases = pdf_font_check(pdf)
print(f"\nPDF font check {pdf}: Type3 = {n3}, ArialMT present = {arial}, BaseFonts = {bases}")
assert n3 == 0 and arial
