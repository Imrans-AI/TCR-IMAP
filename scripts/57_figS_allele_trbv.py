#!/usr/bin/env python3
"""SI figure S8 (p84v2) - per-complex pairwise AUPRC under tcr70, grouped by MHC allele and by TRBV. Descriptive only.

A  by MHC allele (complex_metadata.allele, regex-parsed from the STCRDab compound field; 108 of 174 complexes carry
   none and are shown as one 'not annotated' group, set apart)
B  by TRBV, as recorded in complex_metadata.trbv (subgroup level, e.g. TRBV6)
Each group: every complex for TCR-IMAP and for the within-fold pair-frequency prior (points), the group mean (line),
and the group n. Groups with fewer than 3 complexes are omitted and listed. No test is performed and no marker carries
the fill = significance meaning; this is not a subgroup analysis. Per-complex AUPRC tracks a complex's contact count
and is not a difficulty measure.

Sources (nothing is trained, predicted or re-scored): results/tables/percomplex_tcr70_paired_p84v2.tsv (pairwise tcr70,
TCR-IMAP 'grip' and pair-frequency prior 'prior'; its allele and trbv columns checked against
data/processed/complex_metadata.tsv), checked against pairwise_ours_p84v2.npz and against the tcr70 row of
protocol_main_p84v2.csv (macro-AUPRC 0.3598 and 0.3090).
"""
import sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, "src")
from viz.style import (apply_style, save_fig, print_validation_table, minus, panel_letters, PAGE_WIDTH_IN, MODEL_LABEL,
                       METHOD_COLORS)
from viz.figstats import pdf_font_check

apply_style()
T = "results/tables/"
OUT_STEM = "figures/figS_allele_trbv_p84v2_v1"
NMIN = 3
TXT = "0.30"
RD = dict(float_precision="round_trip")
CM, CP = METHOD_COLORS[MODEL_LABEL], METHOD_COLORS["Pair-frequency prior"]

pr = pd.read_csv(T + "percomplex_tcr70_paired_p84v2.tsv", sep="\t", **RD)
po = np.load(T + "pairwise_ours_p84v2.npz", allow_pickle=True)
assert (po["pdb"] == pr.pdb.values).all() and np.allclose(po["auprc"], pr.grip.values) and len(pr) == 174
md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
for c in ("allele", "trbv"):
    assert (md.loc[pr.pdb, c].fillna("NA").values == pr[c].fillna("NA").values).all(), c
pm = pd.read_csv(T + "protocol_main_p84v2.csv").set_index("protocol").loc["tcr70"]
assert round(pr.grip.mean(), 4) == pm.model_auprc and round(pr.prior.mean(), 4) == pm.freq_auprc
assert pr[["grip", "prior"]].notna().all().all()
print(f"174 complexes; macro-AUPRC {pr.grip.mean():.4f} / {pr.prior.mean():.4f} = protocol_main_p84v2.csv tcr70; "
      "allele and trbv columns equal complex_metadata.tsv")

VAL, OMIT = [], {}
def groups(col, na_label):
    g = pr[col].fillna(na_label)
    cnt = g.value_counts()
    keep = [k for k in cnt.index if cnt[k] >= NMIN and k != na_label]
    keep = sorted(keep, key=lambda k: (-cnt[k], k))
    om = cnt[[k for k in cnt.index if cnt[k] < NMIN]]
    OMIT[col] = om
    assert cnt.sum() == 174
    return g, cnt, keep

FW, FH = PAGE_WIDTH_IN, 5.6
fig = plt.figure(figsize=(FW, FH))
rng = np.random.default_rng(0)        # jitter only
axes = {}
for lab, col, na, (b, h) in (("A", "allele", "not annotated", (3.35, 1.65)), ("B", "trbv", "NA", (0.95, 1.65))):
    g, cnt, keep = groups(col, na)
    order = keep + ([na] if col == "allele" else [])
    ax = fig.add_axes([0.62 / FW, b / FH, 6.45 / FW, h / FH]); axes[lab] = ax
    xs = {}
    for i, k in enumerate(order):
        x = i + (0.6 if k == na else 0.0)
        xs[k] = x
        sel = (g == k).values
        for dx, v, c in ((-0.18, pr.grip.values[sel], CM), (+0.18, pr.prior.values[sel], CP)):
            ax.scatter(x + dx + rng.uniform(-0.09, 0.09, sel.sum()), v, s=4, color=c, alpha=0.55, edgecolor="none", zorder=2)
            ax.plot([x + dx - 0.15, x + dx + 0.15], [v.mean(), v.mean()], color=c, lw=1.6, zorder=3, solid_capstyle="butt")
        ax.text(x, 1.0, f"n = {sel.sum()}", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.8,
                color=TXT)
        VAL.append(dict(panel=lab, group=k, n=int(sel.sum()), m=f"{pr.grip.values[sel].mean():.4f}",
                        p=f"{pr.prior.values[sel].mean():.4f}"))
    if col == "allele":
        ax.axvline((xs[order[-2]] + xs[na]) / 2, color="0.6", lw=0.6, ls=(0, (2, 2)))
    ax.set_xticks([xs[k] for k in order])
    ax.set_xticklabels([k.replace("HLA-", "") for k in order], fontsize=6.3,
                       rotation=0 if col == "allele" else 35, ha="center" if col == "allele" else "right",
                       rotation_mode="anchor")
    ax.set_xlim(-0.6, max(xs.values()) + 0.6); ax.set_ylim(0, 1.0)
    ax.set_ylabel("Pairwise AUPRC\n(per complex, tcr70)", linespacing=1.2)
    om = OMIT[col]
    note = (f"omitted (n < {NMIN}): {len(om)} {'alleles' if col == 'allele' else 'TRBV groups'}, {int(om.sum())} complexes"
            + (" (incl. 1 without TRBV)" if col == "trbv" and "NA" in om.index else ""))
    ax.set_title(f"By MHC allele" if col == "allele" else "By TRBV", fontsize=8, pad=13, loc="left")
    ax.text(1.0, 1.10, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8, color=TXT)
fig.legend(handles=[Line2D([], [], marker="o", ls="", color=CM, ms=3.5, label=f"{MODEL_LABEL}, per complex"),
                    Line2D([], [], color=CM, lw=1.6, label=f"{MODEL_LABEL}, group mean"),
                    Line2D([], [], marker="o", ls="", color=CP, ms=3.5, label="pair-frequency prior, per complex"),
                    Line2D([], [], color=CP, lw=1.6, label="pair-frequency prior, group mean")],
           loc="lower center", bbox_to_anchor=(0.5, 0.005), ncol=4, frameon=False, fontsize=6.3, handletextpad=0.4,
           columnspacing=1.4)
panel_letters(fig, {"A": [axes["A"]], "B": [axes["B"]]})
png, pdf = save_fig(fig, OUT_STEM)
plt.close(fig)

print("\n--- groups drawn (n >= 3), group means; source percomplex_tcr70_paired_p84v2.tsv (grip, prior) ---")
print_validation_table(VAL, [("panel", "panel", None), ("group", "group", None), ("n", "n", "%d"),
                             ("m", "TCR-IMAP mean", None), ("p", "prior mean", None)])
for col, om in OMIT.items():
    print(f"omitted {col} groups (n < {NMIN}): " + ", ".join(f"{k} ({v})" for k, v in om.items()))
n3, arial, bases = pdf_font_check(pdf)
print(f"\nPDF font check {pdf}: Type3 = {n3}, ArialMT present = {arial}, BaseFonts = {bases}")
assert n3 == 0 and arial
