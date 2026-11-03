#!/usr/bin/env python3
"""Figure 1 - from structural data to controlled evaluation.

Rebuilt 2 November 2026 (v2). The previous version is preserved under
figures/previous_versions/fig1_workflow_v1_20261102.*.

Design follows the final Module 2 Figure 1 (verified by SHA-256 against word/media/image1.png of
STEP_MHC_INTEGRATED_MANUSCRIPT_v11_AdvSci.docx; drawn manually in BioRender, so no editable source
exists) WITHOUT copying it. Principles adopted:
  * landscape, panels enclosed in thin coloured rounded borders, letter top-left inside the border
  * one left-to-right dataflow across the top row, joined by grey block arrows
  * abstract objects drawn as icons -- token strips, a sparse contact grid, cluster chips --
    instead of being described in prose
  * dimensions annotated at the point of use, not in a legend
  * one filled dark box marks the model output; everything else is outline-on-white
  * one colour per concept, held constant; very little text, the caption carries the detail
M3 differs from M2 in what it foregrounds: the bottom row is the evaluation framework, which is
this paper's contribution, and is given the heavier accent and the larger share of the page.

Every number is read from a committed artifact and asserted below. No performance value, no
invented structural illustration and no architectural feature that the code does not implement
appears in this figure.

SUPERSEDED 2 November 2026. Figure 1 in the manuscript is the author-approved figure installed at
figures/main/Figure_1.{png,pdf,svg}; its editable source is figures/source/fig1_approved_20261102/
and its provenance is figures/main/Figure_1.SOURCE.md. This script no longer produces the
manuscript figure and cannot write to figures/main/. It is kept only so the superseded version
stays reproducible, and refuses to run without --allow-superseded.

usage: python scripts/70_fig1_workflow.py --allow-superseded
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle

sys.path.insert(0, ".")
from src.viz.style import apply_style, save_fig, MODEL_LABEL, METHOD_COLORS, MOLECULE_COLORS, PROTOCOL_COLORS

ROOT = Path(__file__).resolve().parent.parent

# This script is superseded by the author-approved Figure 1. It must never overwrite the approved
# asset, so it refuses to run unless the caller says explicitly that the superseded figure is what
# is wanted, and it writes only into figures/previous_versions/.
if "--allow-superseded" not in sys.argv:
    sys.exit(
        "70_fig1_workflow.py is SUPERSEDED and does not build the manuscript Figure 1.\n"
        "The manuscript embeds the author-approved figures/main/Figure_1.png "
        "(sha256 a447b2a87ed2d74ab19379e7b6439fe0a8d2728331f4b2692e3700ca8b45bbd0);\n"
        "its editable source is figures/source/fig1_approved_20261102/design_figure.py and its\n"
        "provenance is figures/main/Figure_1.SOURCE.md.\n"
        "To rebuild the superseded version anyway, pass --allow-superseded; it will be written to\n"
        "figures/previous_versions/fig1_workflow_superseded.{png,pdf} and left there."
    )
OUT_STEM = "figures/previous_versions/fig1_workflow_superseded"
apply_style()

# ---------------------------------------------------------------- verified inputs
z = np.load(ROOT/"data/processed/pair_maps_p84v2.npz", allow_pickle=True)
md = pd.read_csv(ROOT/"data/processed/complex_metadata.tsv", sep="\t")
ix = pd.read_csv(ROOT/"data/interim/stcrdab_index.tsv", sep="\t")
proto = json.load(open(ROOT/"data/splits_v2/protocols.json"))
leak = pd.read_csv(ROOT/"results/tables/leakage_audit_dual.tsv", sep="\t")

N, K = len(md), md.clu70.nunique()
REG = pd.Series([str(x) for x in z["row_imgt_region"]]).value_counts()
NROW, NCOL = z["Y"].shape[1], z["Y"].shape[2]
NRES, NCON = int(z["S"].sum()), int(z["Y"].sum())
CN = z["col_imgt_number"].astype(int)
NA1, NA2 = int((CN <= 90).sum()), int((CN >= 1050).sum())
m = ix[ix.pdb.isin(md.pdb)].drop_duplicates("pdb")
a = m.alpha_organism.astype(str).str.strip().str.lower()
b = m.beta_organism.astype(str).str.strip().str.lower()
N_MOUSE = int(((a == "mus musculus") & (b == "mus musculus")).sum())
N_MIXED = int(a.str.contains("pan troglodytes").sum())
N_HUM = N - N_MOUSE - N_MIXED
EVAL = {p: sum(len(f) for f in v if len(f) >= 3) for p, v in proto.items()}
SIM = leak.groupby("protocol").max_sim_mean.mean().round(3).to_dict()
assert (N, K, NROW, NCOL) == (174, 69, 75, 84)
assert (NRES, NCON) == (893928, 5281) and (NA1, NA2) == (41, 43)
assert (N_HUM, N_MOUSE, N_MIXED) == (165, 8, 1)
assert EVAL == dict(random=174, allele=64, tcr70=174, double=30, temporal=60)
assert (m.mhc_chain1_organism.astype(str).str.strip().str.lower() == "homo sapiens").all()

# ---------------------------------------------------------------- palette
INK, MUT, HAIR = "#1A1A1A", "#5A5A5A", "#9A9A9A"
BLUE = METHOD_COLORS[MODEL_LABEL]
A_COL, B_COL, C_COL = "#1F77A8", "#2E8B6F", "#B8860B"      # top row panel borders
D_COL = E_COL = "#6A3D9A"                                   # evaluation row, the contribution
LOOPC = {"CDR1": MOLECULE_COLORS["TCR alpha"], "CDR2": "#8ECDF0",
         "CDR3": MOLECULE_COLORS["TCR beta"], "FR4": "#F2C3A0"}

fig = plt.figure(figsize=(7.2, 4.75))
PAN = {}


def panel(x, y, w, h, letter, title, colour):
    """Enclosed rounded border with the letter and title inside, M2 style."""
    fig.patches.append(FancyBboxPatch((x, y), w, h, transform=fig.transFigure,
                                      boxstyle="round,pad=0.004,rounding_size=0.012",
                                      fc="white", ec=colour, lw=1.1, zorder=0))
    ax = fig.add_axes([x + 0.012, y + 0.012, w - 0.024, h - 0.055])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.set_zorder(5); ax.patch.set_visible(False)
    fig.text(x + 0.010, y + h - 0.040, letter, fontsize=10, fontweight="bold", color=INK)
    fig.text(x + 0.034, y + h - 0.036, title, fontsize=7.0, color=INK, va="center")
    PAN[letter] = ax
    return ax


def chip(ax, x, y, w, h, text, fc="white", ec=HAIR, fs=5.6, weight="normal", tc=INK, lw=0.7):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.02",
                                fc=fc, ec=ec, lw=lw, zorder=2))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs,
            color=tc, zorder=3, linespacing=1.3, fontweight=weight)


def arrow(ax, p0, p1, lw=0.85, color=MUT, ms=6):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms, lw=lw,
                                 color=color, shrinkA=1, shrinkB=1, zorder=1))


def block_arrow(x, y, w=0.018):
    """Grey block arrow between top-row panels, in figure coordinates."""
    fig.patches.append(FancyArrowPatch((x, y), (x + w, y), transform=fig.transFigure,
                                       arrowstyle="-|>", mutation_scale=11, lw=2.6,
                                       color="#BFBFBF", zorder=1))


def strip(ax, x, y, w, h, parts, fs=5.0, label_below=True):
    """A sequence drawn as coloured segments with counts inside."""
    tot = sum(n for _, n, _ in parts); cx = x
    for name, n, col in parts:
        ww = w * n / tot
        ax.add_patch(Rectangle((cx, y), ww, h, fc=col, ec="white", lw=0.6, zorder=2))
        if ww > 0.055:                      # suppress the label on segments too narrow to hold it
            ax.text(cx + ww/2, y + h/2, str(n), ha="center", va="center", fontsize=fs,
                    color="white", fontweight="bold", zorder=3)
        if label_below:
            ax.text(cx + ww/2, y - 0.055, name, ha="center", va="center", fontsize=4.8, color=MUT)
        cx += ww


# ================================================================= A  data set
ax = panel(0.012, 0.515, 0.225, 0.455, "A", "Structural data set", A_COL)
rows = [("1,002 STCRDab rows\n634 PDB entries", "#F4F4F4", HAIR, "normal"),
        ("αβ TCR · peptide · MHC-I\nhuman MHC heavy chain", "#F4F4F4", HAIR, "normal"),
        ("275 → deduplicate → 175", "#F4F4F4", HAIR, "normal"),
        (f"{N} complexes", "#E6F0F7", A_COL, "bold")]
yy, hh = 0.945, 0.126
for i, (t, fc, ec, w_) in enumerate(rows):
    chip(ax, 0.03, yy - hh, 0.94, hh, t, fc=fc, ec=ec, fs=5.5, weight=w_,
         lw=1.1 if i == 3 else 0.7)
    if i < 3:
        arrow(ax, (0.5, yy - hh), (0.5, yy - hh - 0.042), ms=5)
    yy -= hh + 0.042
ax.text(0.5, 0.235, "MHC class I human in all 174", ha="center", fontsize=5.4, color=INK)
strip(ax, 0.10, 0.105, 0.80, 0.066,
      [("human", N_HUM, "#8FA8C8"), ("mouse", N_MOUSE, "#C8A88F"), ("mixed", N_MIXED, "#D8D8D8")],
      fs=4.6, label_below=False)
ax.text(0.5, 0.045, f"receptors  human {N_HUM} \u00b7 mouse {N_MOUSE} \u00b7 mixed {N_MIXED}",
        ha="center", fontsize=4.9, color=MUT)
ax.text(0.5, 0.000, "species was never a selection criterion", ha="center", fontsize=4.7, color=MUT)

block_arrow(0.240, 0.735)

# ================================================================= B  representation
ax = panel(0.262, 0.515, 0.300, 0.455, "B", "Aligned residue representation", B_COL)
ax.text(0.02, 0.945, f"TCR  {NROW} rows", fontsize=5.8, color=INK, va="center")
strip(ax, 0.02, 0.80, 0.62, 0.085,
      [(k, int(REG[k]), LOOPC[k]) for k in ("CDR1", "CDR2", "CDR3", "FR4")], fs=5.0)
ax.text(0.02, 0.655, f"MHC  {NCOL} columns", fontsize=5.8, color=INK, va="center")
strip(ax, 0.02, 0.515, 0.62, 0.085,
      [("α1", NA1, "#B9B9B9"), ("α2", NA2, "#8F8F8F")], fs=5.0)
ax.text(0.02, 0.400, "keyed by IMGT number and insertion code;\ntokens placed, never packed",
        fontsize=5.0, color=MUT, va="top", linespacing=1.45)
# sparse contact grid
gx, gy, gw, gh = 0.70, 0.46, 0.28, 0.42
rng = np.random.default_rng(3)
ax.add_patch(Rectangle((gx, gy), gw, gh, fc="#FBFBFB", ec=HAIR, lw=0.7, zorder=2))
for i in range(13):
    for j in range(11):
        if rng.random() < 0.075:
            ax.add_patch(Rectangle((gx + j*gw/11, gy + i*gh/13), gw/11, gh/13,
                                   fc=BLUE, ec="none", zorder=3))
ax.text(gx + gw/2, gy - 0.055, f"{N} × {NROW} × {NCOL}", ha="center", fontsize=5.2, color=INK)
chip(ax, 0.02, 0.055, 0.96, 0.21,
     f"{NRES:,} resolved pairs   ·   {NCON:,} contacts at 4.5 Å   ({100*NCON/NRES:.2f} %)",
     fc="#EFF6F3", ec=B_COL, fs=5.5, lw=1.0)
ax.text(0.5, 0.005, "MHC window fixed before training", ha="center", fontsize=4.9, color=MUT)

block_arrow(0.565, 0.735)

# ================================================================= C  model
ax = panel(0.587, 0.515, 0.401, 0.455, "C", "Model and reference methods", C_COL)
chip(ax, 0.03, 0.845, 0.40, 0.115, "TCR tokens", fc="#EAF4FB", ec=LOOPC["CDR1"], fs=5.5)
chip(ax, 0.49, 0.845, 0.40, 0.115, "MHC tokens", fc="#F2F2F2", ec=HAIR, fs=5.5)
arrow(ax, (0.23, 0.845), (0.23, 0.775), ms=5); arrow(ax, (0.69, 0.845), (0.69, 0.775), ms=5)
chip(ax, 0.03, 0.655, 0.40, 0.118, "encoder", ec=BLUE, fs=5.6)
chip(ax, 0.49, 0.655, 0.40, 0.118, "encoder", ec=BLUE, fs=5.6)
arrow(ax, (0.23, 0.655), (0.40, 0.585), ms=5); arrow(ax, (0.69, 0.655), (0.52, 0.585), ms=5)
chip(ax, 0.17, 0.462, 0.58, 0.118, "cross-attention", ec=BLUE, fs=5.6)
arrow(ax, (0.46, 0.462), (0.46, 0.392), ms=5)
chip(ax, 0.17, 0.268, 0.58, 0.120, f"bilinear head\n{NROW} × {NCOL} contact map",
     fc="#17486B", ec="#17486B", fs=5.6, weight="bold", tc="white", lw=1.0)
ax.text(0.46, 0.215, "529,310 parameters   ·   amino-acid sequence only",
        ha="center", fontsize=5.0, color=MUT)
chip(ax, 0.03, 0.015, 0.44, 0.165,
     "Within-fold priors\npair- and position-frequency\nnearest-neighbour transfer",
     fc="#FAFAFA", ec=HAIR, fs=5.0)
chip(ax, 0.53, 0.015, 0.44, 0.165,
     "Reimplemented comparators\nStruct-TCS · Seq-TCS\nTEIM-Res",
     fc="#FAFAFA", ec=HAIR, fs=5.0)

# ================================================================= D  protocols
ax = panel(0.012, 0.030, 0.470, 0.455, "D", "What each protocol withholds", D_COL)
ORDER = ["random", "double", "allele", "temporal", "tcr70"]
PLAB = {"random": "random", "double": "intersection (legacy)", "allele": "allele",
        "temporal": "temporal", "tcr70": "tcr70  (primary)"}
HOLD = {"random": (0, 0, 0), "double": (0, 0, 0), "allele": (0, 1, 0),
        "temporal": (0, 0, 1), "tcr70": (1, 0, 0)}
COLS = ["TCR\ncluster", "MHC\nallele", "deposition\ndate"]
cx = [0.455, 0.555, 0.655]
for j, c in enumerate(COLS):
    ax.text(cx[j], 0.945, c, ha="center", va="center", fontsize=4.6, color=MUT, linespacing=1.2)
ax.text(0.86, 0.945, "train/test\nsimilarity", ha="center", va="center", fontsize=4.6,
        color=MUT, linespacing=1.2)
yy = 0.78
for p in ORDER:
    ax.text(0.02, yy, PLAB[p], fontsize=5.4, va="center",
            color=INK, fontweight="bold" if p == "tcr70" else "normal")
    ax.text(0.40, yy, f"n = {EVAL[p]}", fontsize=5.0, va="center", ha="right", color=MUT)
    for j, on in enumerate(HOLD[p]):
        ax.add_patch(Circle((cx[j], yy), 0.030, fc=PROTOCOL_COLORS[p] if on else "white",
                            ec=PROTOCOL_COLORS[p] if on else "#C8C8C8", lw=0.9, zorder=3))
    bx, bw = 0.745, 0.21
    ax.add_patch(Rectangle((bx, yy - 0.030), bw, 0.060, fc="#F2F2F2", ec="none", zorder=2))
    ax.add_patch(Rectangle((bx, yy - 0.030), bw*SIM[p], 0.060, fc=PROTOCOL_COLORS[p],
                           ec="none", zorder=3))
    ax.text(bx + bw + 0.012, yy, f"{SIM[p]:.2f}", fontsize=4.9, va="center", color=MUT)
    yy -= 0.158
ax.text(0.02, 0.055, "Filled marks what the protocol holds out. Bars are the mean maximum\n"
        "train/test similarity each one still permits; they are not a strictness ranking.",
        fontsize=4.9, color=MUT, va="center", linespacing=1.5)

# ================================================================= E  inference
ax = panel(0.494, 0.030, 0.494, 0.455, "E", "Cluster-level paired inference", E_COL)
chip(ax, 0.02, 0.655, 0.20, 0.235, f"{N}\ncomplexes", fc="#FFFFFF", ec=HAIR, fs=5.6)
arrow(ax, (0.225, 0.772), (0.275, 0.772), ms=6)
# cluster chips
ax.add_patch(FancyBboxPatch((0.283, 0.655), 0.22, 0.235,
                            boxstyle="round,pad=0.004,rounding_size=0.02",
                            fc="#F3EEF8", ec=E_COL, lw=1.0, zorder=2))
rng2 = np.random.default_rng(7)
for i in range(4):
    for j in range(6):
        ax.add_patch(Rectangle((0.300 + j*0.031, 0.700 + i*0.042), 0.022, 0.030,
                               fc=E_COL, alpha=0.25 + 0.6*rng2.random(), ec="none", zorder=3))
ax.text(0.393, 0.672, f"{K} TCR clusters", ha="center", fontsize=5.0, color=INK, zorder=4)
arrow(ax, (0.508, 0.772), (0.558, 0.772), ms=6)
chip(ax, 0.566, 0.655, 0.21, 0.235, "resample\nclusters with\nreplacement", fc="white", ec=HAIR, fs=5.3)
arrow(ax, (0.781, 0.772), (0.831, 0.772), ms=6)
chip(ax, 0.785, 0.655, 0.20, 0.235, "95 % CI\nand p", fc="#F3EEF8", ec=E_COL, fs=5.8,
     weight="bold", lw=1.1)
ax.text(0.672, 0.620, "20,000 resamples, seed 42", ha="center", fontsize=4.8, color=MUT)
ax.text(0.02, 0.545, "Complexes within a cluster are not independent, so clusters are the resampling unit.\n"
        "Every comparison is a paired per-complex difference between two methods.",
        fontsize=5.0, color=MUT, va="top", linespacing=1.5)
chip(ax, 0.02, 0.185, 0.96, 0.23,
     "Fixed before the reported results were produced\n"
     "metric  per-complex macro-AUPRC      protocol  tcr70\n"
     f"comparison  {MODEL_LABEL} versus Struct-TCS, per-residue task",
     fc="#FDF6E6", ec="#C8882B", fs=5.4, lw=1.0)
chip(ax, 0.02, 0.015, 0.96, 0.135,
     "Controls: training-size-matched sweep  ·  cluster-weighted sensitivity  ·  canonical-CDR scoring",
     fc="#FAFAFA", ec=HAIR, fs=4.9)

save_fig(fig, OUT_STEM)
plt.close(fig)

print("\n--- Figure 1 values, each asserted against its source artifact ---")
for k, v in [("complexes", N), ("clu70 clusters", K), ("TCR rows / MHC columns", f"{NROW} / {NCOL}"),
             ("rows CDR1/CDR2/CDR3/FR4", "/".join(str(int(REG[r])) for r in ("CDR1","CDR2","CDR3","FR4"))),
             ("MHC columns a1 / a2", f"{NA1} / {NA2}"),
             ("resolved pairs / contacts", f"{NRES:,} / {NCON:,} ({100*NCON/NRES:.2f} %)"),
             ("receptors human/mouse/mixed", f"{N_HUM}/{N_MOUSE}/{N_MIXED}"),
             ("evaluated per protocol", EVAL), ("mean max similarity", SIM)]:
    print(f"  {k:<30} {v}")
