#!/usr/bin/env python3
"""Figure 2 (p82v1) - dataset redundancy and leakage.

Port of scripts/38_fig2_dataset.py onto src/viz/style.py (apply_style, Type 42 fonts, paper
palette). Inputs are TCR-sequence and split based, so they are unaffected by the 24 Aug
candidate-set correction; no MHC position enters this figure.

A  V-gene usage, top 8 TRAV / TRBV (colour = TCR chain, MOLECULE_COLORS)
B  all-pairs vs nearest-neighbour TCR Levenshtein similarity ratio (174 complexes)
C  complexes by cluster size at 70 / 80 / 90 % (sequential CMAP_SEQ shades, darker = stricter)
D  per-complex max similarity to any training complex, per protocol (PROTOCOL_COLORS),
   curves labelled directly at their right-hand ends (no legend; line style is reserved for loop class)
E  split duplicate families per fold, allele-independent definition (CDR3 pair + peptide)

Allele coverage fix (HANDOFF_1.md section 6, known bug 1): 66 complexes are assigned an allele
fold but 64 are evaluated -- fold 3 (n = 2) is dropped by the frozen pipeline's n < 3 guard.
That fold is dropped here too, in both D and E, so every allele number is n = 64, matching
the evaluation tables. The script asserts that the remaining 64 equal the evaluated set in
percomplex_allele_p82v1.npz.

Similarity is Levenshtein.ratio over concatenated CDR1/2/3 of both chains, an edit-distance
similarity -- labelled "Levenshtein similarity ratio", never identity.
Sources: data/processed/complex_metadata.tsv, results/tables/tcr_identity_matrix.npy,
results/tables/percomplex_maxsim.tsv, results/tables/leakage_audit_dual.tsv,
results/tables/percomplex_{protocol}_p82v1.npz (evaluated-set check only).
"""
import sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
sys.path.insert(0, "src")
from viz.style import (apply_style, save_fig, print_validation_table, PAGE_WIDTH_IN,
                       PANEL_LABEL_STYLE as PL, panel_letters, SIM_LABEL, MOLECULE_COLORS, PROTOCOL_COLORS,
                       CMAP_SEQ, minus)
from viz.figstats import pdf_font_check

apply_style()
PANELS = {}

OUT_STEM = "figures/fig2_dataset_p82v1_v1"
ORDER = ["random", "double", "allele", "temporal", "tcr70"]   # keys as stored in the frozen files
# Display labels. "double" is a legacy intersection-selected partition, not a double holdout:
# its held-out groups remain in training (results/review/V3_PART_A_VERIFICATION.md A1), so it is
# named accordingly and the protocols are NOT presented as an ordering by strictness.
PLAB = {"random": "random", "double": "intersection\n(legacy)", "allele": "allele",
        "temporal": "temporal", "tcr70": "tcr70"}
C_A, C_B = MOLECULE_COLORS["TCR alpha"], MOLECULE_COLORS["TCR beta"]
C_ALLPAIRS, C_NN = CMAP_SEQ(0.40), CMAP_SEQ(1.0)
C_CLU = {"clu70": CMAP_SEQ(1.0), "clu80": CMAP_SEQ(0.62), "clu90": CMAP_SEQ(0.36)}
TXT = "0.30"   # neutral annotation text, as in the Fig 3/5 templates

m = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t")
pc = pd.read_csv("results/tables/percomplex_maxsim.tsv", sep="\t")
au = pd.read_csv("results/tables/leakage_audit_dual.tsv", sep="\t")
S = np.load("results/tables/tcr_identity_matrix.npy")
assert len(m) == 174 and S.shape == (174, 174) and "6uz1" not in set(m.pdb)

# ---- allele coverage fix: drop the allele fold the frozen pipeline never evaluates -------
fold_n = pc[pc.protocol == "allele"].groupby("fold").size()
small = fold_n[fold_n < 3].index.tolist()
assert small == [3] and int(fold_n[3]) == 2, fold_n.to_dict()
dropped = pc[(pc.protocol == "allele") & pc.fold.isin(small)].pdb.tolist()
pc = pc[~((pc.protocol == "allele") & pc.fold.isin(small))].copy()
au = au[~((au.protocol == "allele") & au.fold.isin(small))].copy()
for p in ORDER:   # every curve in D is exactly the evaluated set of that protocol
    ev = set(np.load(f"results/tables/percomplex_{p}_p82v1.npz", allow_pickle=True)["pdb"])
    assert set(pc.loc[pc.protocol == p, "pdb"]) == ev, p
    assert int(au.loc[au.protocol == p, "n_test"].sum()) == len(ev), p
N_EVAL = pc.groupby("protocol").size().to_dict()
assert N_EVAL["allele"] == 64

iu = np.triu_indices_from(S, k=1); off = S[iu]
D = S.copy(); np.fill_diagonal(D, -1); nn = D.max(axis=1)

fig = plt.figure(figsize=(PAGE_WIDTH_IN, 5.0), constrained_layout=True)
gs = fig.add_gridspec(2, 6)

# --- A  V-gene usage --------------------------------------------------------------------
ax = fig.add_subplot(gs[0, 0:2])
ta, tb = m.trav.value_counts().head(8), m.trbv.value_counts().head(8)
labs = list(ta.index) + [""] + list(tb.index)
vals = list(ta.values) + [0] + list(tb.values)
cols = [C_A] * 8 + ["none"] + [C_B] * 8
y = np.arange(len(labs))[::-1]
ax.barh(y, vals, color=cols, height=0.72, linewidth=0)
ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=6.2)
ax.set_xlabel("Complexes"); ax.set_xlim(0, 72)
for yy, v in zip(y, vals):
    if v: ax.text(v + 1.5, yy, minus(str(v)), va="center", fontsize=5.8, color=TXT)
ax.text(44, 12.6, "TRAV", color=C_A, fontsize=7.5, fontweight="bold", va="center")
ax.text(44, 3.4, "TRBV", color=C_B, fontsize=7.5, fontweight="bold", va="center")
ax.set_title("V-gene usage (top 8)", fontsize=8, pad=3)
PANELS["A"] = [ax]

# --- B  pairwise vs nearest-neighbour similarity ------------------------------------------
ax = fig.add_subplot(gs[0, 2:4])
for v, c, lw, lab in ((off, C_ALLPAIRS, 1.1, "All pairs (n = %s)" % f"{len(off):,}"),
                      (nn, C_NN, 1.6, "Nearest neighbor (n = %d)" % len(nn))):
    sv = np.sort(v)
    ax.step(sv, np.arange(1, len(sv) + 1) / len(sv), where="post", color=c, lw=lw, label=lab)
ax.axvline(0.90, ls="--", lw=0.8, color="0.4")
ax.set_xlabel(f"TCR {SIM_LABEL}"); ax.set_ylabel("Cumulative fraction")
ax.set_xlim(0.25, 1.02); ax.set_ylim(0, 1.03)
ax.legend(frameon=False, fontsize=5.8, loc="upper left", handlelength=1.1)
n_exact = int((nn >= 0.999).sum()); frac_exact = n_exact / len(nn)
# set vertically between the dashed 0.90 line and the step at 1.0 it describes, so no
# leader line has to cross the dashed line
ax.text(0.951, 0.72, minus("%.1f%% (%d/174)\nhave an exact twin" % (100 * frac_exact, n_exact)),
        rotation=90, fontsize=6, color=C_NN, ha="center", va="center", linespacing=1.3)
ax.text(0.03, 0.80, minus("ratio ≥ 0.90:\n%.1f%% of pairs\n%.1f%% of complexes"
        % (100 * (off >= 0.9).mean(), 100 * (nn >= 0.9).mean())),
        transform=ax.transAxes, fontsize=6, va="top", ha="left", linespacing=1.45, color=TXT)
ax.set_title("Redundancy: pairs vs neighbors", fontsize=8, pad=3)
PANELS["B"] = [ax]

# --- C  cluster-size mass -----------------------------------------------------------------
ax = fig.add_subplot(gs[0, 4:6])
BL = ["1", "2", "3", "4-5", "6-10", ">10"]
def b(s): return 0 if s == 1 else 1 if s == 2 else 2 if s == 3 else 3 if s <= 5 else 4 if s <= 10 else 5
w = 0.26
n_clu = {}
for k, col in enumerate(["clu70", "clu80", "clu90"]):
    vc = m[col].value_counts(); mass = np.zeros(6); n_clu[col] = len(vc)
    for _, sz in vc.items(): mass[b(sz)] += sz
    ax.bar(np.arange(6) + (k - 1) * w, mass, w, color=C_CLU[col], linewidth=0,
           label="%d%% (%d clusters)" % ({"clu70": 70, "clu80": 80, "clu90": 90}[col], len(vc)))
ax.set_xticks(range(6)); ax.set_xticklabels(BL, fontsize=6.5)
ax.set_xlabel("Cluster size"); ax.set_ylabel("Complexes"); ax.set_ylim(0, 72)
ax.legend(frameon=False, fontsize=5.8, loc="upper center", ncol=1, handlelength=0.9, borderpad=0.2)
ax.set_title("Complexes by cluster size", fontsize=8, pad=3)
PANELS["C"] = [ax]

# --- D  per-complex max similarity to training --------------------------------------------
ax = fig.add_subplot(gs[1, 0:3])
end_y = {}
for p in ORDER:
    v = np.sort(pc.loc[pc.protocol == p, "max_sim"].values)
    ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=PROTOCOL_COLORS[p],
            lw=1.7 if p in ("tcr70", "random") else 1.2)
    end_y[p] = float((v < 0.9995).mean())   # plateau just left of the final step at 1.0
ax.axvline(0.90, ls="--", lw=0.8, color="0.4")
ax.set_xlabel(f"Max {SIM_LABEL} to any training complex")
ax.set_ylabel("Cumulative fraction"); ax.set_xlim(0.52, 1.14); ax.set_ylim(0, 1.03)
ax.set_xticks([0.6, 0.7, 0.8, 0.9, 1.0]); ax.spines["bottom"].set_bounds(0.52, 1.0)
r_ms = pc.loc[pc.protocol == "random", "max_sim"]; t_ms = pc.loc[pc.protocol == "tcr70", "max_sim"]
# direct labels (no legend, no line styles): the four curves that end at 1.0 are labelled just
# right of their final step at their plateau height, pushed apart where they crowd; tcr70 ends
# at its maximum and is labelled to the left of its top end
ends = sorted((end_y[p], p) for p in ORDER if p != "tcr70")
placed, GAPY = [], 0.085
for yv, p in ends:
    yl = max(yv, placed[-1][0] + GAPY) if placed else yv
    placed.append((yl, p))
over = placed[-1][0] - 0.99
if over > 0: placed = [(yl - over, p) for yl, p in placed]
LABEL_Y = {p: yl for yl, p in placed}
for p in ORDER:
    if p == "tcr70":
        ax.text(t_ms.max() - 0.007, 0.985, f"tcr70 (n = {N_EVAL[p]})", ha="right", va="center",
                fontsize=6, color=PROTOCOL_COLORS[p])
    else:
        ax.text(1.008, LABEL_Y[p], f"{PLAB[p].replace(chr(10), ' ')} (n = {N_EVAL[p]})",
                ha="left", va="center", fontsize=6, color=PROTOCOL_COLORS[p])
ax.text(0.886, 0.03, minus("random: median %.3f,\n%.0f%% ≥ 0.90\ntcr70: max %.3f,\n%.0f%% ≥ 0.90"
        % (np.median(r_ms), 100 * (r_ms >= 0.9).mean(), t_ms.max(), 100 * (t_ms >= 0.9).mean())),
        transform=blended_transform_factory(ax.transData, ax.transAxes), fontsize=6,
        va="bottom", ha="right", linespacing=1.3, color=TXT)
ax.set_title("Test-set proximity to training", fontsize=8, pad=3)
PANELS["D"] = [ax]

# --- E  duplicate-family leakage ------------------------------------------------------------
ax = fig.add_subplot(gs[1, 3:6])
g = au.groupby("protocol").agg(n=("n_test", "sum"), folds=("fold", "nunique"),
                               lc=("leak_conservative", "mean")).reindex(ORDER)
x = np.arange(len(ORDER))
ax.bar(x, g.lc, 0.55, color=[PROTOCOL_COLORS[p] for p in ORDER])
for xi, v in zip(x, g.lc):
    ax.text(xi, v + 0.16, minus(f"{v:.2f}" if v % 1 and round(v, 1) != v else f"{v:.1f}"),
            ha="center", fontsize=6.5, color="0.25")
ax.set_xticks(x)
ax.set_xticklabels(["%s\n%d/174" % (PLAB[p], g.n[p]) for p in ORDER], fontsize=5.8, linespacing=1.3)
ax.set_ylabel("Split duplicate families per fold")
ax.set_ylim(0, max(g.lc) * 1.25)
ax.set_title("Duplicate-family leakage", fontsize=8, pad=3)
ax.text(0.02, 0.90, "complexes evaluated shown below each protocol",
        transform=ax.transAxes, fontsize=5.8, color="0.45")
PANELS["E"] = [ax]

panel_letters(fig, PANELS)          # one fixed offset for every figure (style.py)
png, pdf = save_fig(fig, OUT_STEM)
plt.close(fig)

# ---------------- panel-to-source validation ----------------------------------------------
vc_ = m.dup_family.value_counts(); ind_ = (m.dup_family.map(vc_) > 1).values
print(f"\nallele fold dropped: fold {small[0]} (n = 2: {', '.join(dropped)}); "
      f"allele n = {N_EVAL['allele']} in D and E, equal to percomplex_allele_p82v1.npz")
rows = [
    dict(panel="B", quantity="complexes", source="complex_metadata.tsv", shown=f"{len(m)}", claude="174"),
    dict(panel="C", quantity="clusters 70/80/90 %", source="complex_metadata.tsv",
         shown=f"{n_clu['clu70']}/{n_clu['clu80']}/{n_clu['clu90']}", claude="69 (70 %); HANDOFF 82/89"),
    dict(panel="B", quantity="exact TCR twin (NN >= 0.999)", source="tcr_identity_matrix.npy",
         shown=f"{n_exact}/174 = {100*frac_exact:.1f}%", claude="57.5 % (100/174)"),
    dict(panel="B", quantity="pairs >= 0.90", source="tcr_identity_matrix.npy",
         shown=f"{100*(off>=0.9).mean():.1f}%", claude="-"),
    dict(panel="B", quantity="complexes NN >= 0.90", source="tcr_identity_matrix.npy",
         shown=f"{100*(nn>=0.9).mean():.1f}%", claude="-"),
    dict(panel="D", quantity="random median max-sim", source="percomplex_maxsim.tsv",
         shown=f"{np.median(r_ms):.3f}", claude="0.983"),
    dict(panel="D", quantity="random >= 0.90", source="percomplex_maxsim.tsv",
         shown=f"{100*(r_ms>=0.9).mean():.1f}%", claude="HANDOFF 60.3 % (fold mean)"),
    dict(panel="D", quantity="tcr70 max max-sim", source="percomplex_maxsim.tsv",
         shown=f"{t_ms.max():.3f}", claude="0.696"),
]
for p in ORDER:
    rows.append(dict(panel="D/E", quantity=f"{p} n evaluated", source="percomplex_maxsim.tsv + leakage_audit_dual.tsv",
                     shown=f"{N_EVAL[p]}", claude={"random": "174", "double": "30", "allele": "64",
                                                    "temporal": "60", "tcr70": "174"}[p]))
for p in ORDER:
    rows.append(dict(panel="E", quantity=f"{p} split dup families/fold ({int(g.folds[p])} folds)",
                     source="leakage_audit_dual.tsv (leak_conservative)", shown=f"{g.lc[p]:.2f}",
                     claude={"random": "7.4", "double": "0.67", "allele": "1.2 (HANDOFF, 5 folds incl. n=2)",
                             "temporal": "0", "tcr70": "0 (no leakage)"}[p]))
print_validation_table(rows, [("panel", "panel", None), ("quantity", "quantity", None),
                              ("source", "source", None), ("shown", "displayed", None),
                              ("claude", "CLAUDE.md / HANDOFF", None)])

# CLAUDE.md "Primary p82v1 results" -> Redundancy line: hard cross-checks
checks = {"complexes == 174": len(m) == 174,
          "clu70 clusters == 69": n_clu["clu70"] == 69,
          "exact twin 57.5 %": round(100 * frac_exact, 1) == 57.5,
          "random median max-sim 0.983": round(float(np.median(r_ms)), 3) == 0.983,
          "tcr70 max max-sim 0.696": round(float(t_ms.max()), 3) == 0.696,
          "no dup-family leakage under tcr70": g.lc["tcr70"] == 0}
print("\nCLAUDE.md cross-check:", {k: ("OK" if v else "MISMATCH") for k, v in checks.items()})
print("\nper-protocol max-sim (pooled per complex):")
print(pc.groupby("protocol").max_sim.agg(["count", "mean", "median", "max",
      lambda s: round((s >= 0.9).mean(), 3)]).rename(columns={"<lambda_0>": "frac>=0.90"})
      .reindex(ORDER).round(3).to_string())
print("median NN: in dup family %.3f (n=%d) | singleton %.3f (n=%d)"
      % (np.median(nn[ind_]), ind_.sum(), np.median(nn[~ind_]), (~ind_).sum()))
print("\npanel D direct labels (curve plateau left of the final step at 1.0 -> label y):")
for p in ORDER:
    if p == "tcr70":
        print(f"  tcr70    ends at max {t_ms.max():.3f}, label left of its top end")
    else:
        print(f"  {p:8s} plateau {end_y[p]:.3f} -> label {LABEL_Y[p]:.3f}")
n3, arial, bases = pdf_font_check(pdf)
print(f"\nPDF font check {pdf}: Type3 = {n3}, ArialMT present = {arial}, BaseFonts = {bases}")
assert n3 == 0 and arial
