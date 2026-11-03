"""Figure 5 -- TCR contact footprint on the MHC alpha1/alpha2 platform (p84v2).

Built from 40c_fig5_geometry_plot_p82v1.py; inputs are the p84v2 descriptive tables, keyed by
(IMGT number, insertion code): 75 TCR rows (68 occupied IMGT positions + 7 CDR3 insertion rows) and
182 G-domain columns, so the window has 84 columns (82 + 1061A + 1072A; alpha2 shows 43).

Rows seen in few complexes. Three CDR3b insertion rows are resolved in fewer than
MIN_ROW_COMPLEXES = 5 complexes (111B and 112C in one, 112B in three). Their frequencies are 1/1
or k/3, not population estimates, and with them the hottest pair would be CDR3b 111B x 1058 at
1.000. They are hatched in A, excluded from A's colour scale (the hottest pair), from the panel B
sum and from the hot-position list (32b-equivalent: scripts/55_pymol_fig5_p84v2.py). D and E count
their contacts as observed (they are real contacts; D and E are rates over resolved pairs).

A  contact frequency, 75 TCR rows x the 84-column window (41 alpha1, 43 alpha2). Strip above
   marks DSSP helix vs non-helical boundary columns. Boxes mark the helix each loop prefers,
   from germline_diagonal_p84v2.csv (cluster bootstrap over 69 clu70 clusters, Holm-corrected).
   v4 drew these from pooled means; they are now cluster-level.
B  MHC platform coloured by contact frequency (PyMOL, 32b). RELATIVE scale: per-position
   frequency summed over all 84 TCR positions, normalised to its own maximum -- NOT panel A's
   per-pair scale.
C  TCR-pMHC complex with CDR loops (PyMOL, 32b)
D  contact rate by loop x MHC region, stating the beta-sheet floor null
E  per-cluster alpha1 - alpha2 contact-rate difference with bootstrap 95% CI

Encoding, per the paper-wide scheme: colour = TCR chain; line style = loop class (germline
solid, CDR3 dashed); fill = significance (filled when the CI excludes zero).

Contact frequency = contacts / resolved complexes per position pair, 4.5 A minimum heavy atom.
Descriptive: no model, no split, no outcome-informed selection. alpha3 (IMGT >= 1091) excluded.
"""
import sys, glob, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import to_rgb, PowerNorm
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
sys.path.insert(0, "src")
from viz.style import (apply_style, save_fig, pstr, PAGE_WIDTH_IN, PANEL_LABEL_STYLE, panel_letters, PANEL_LETTER_OFFSET_IN, PANEL_LETTER_HEIGHT_IN,
                       CMAP_SEQ, REGION_COLORS, MOLECULE_COLORS, GERMLINE_LOOPS,
                       loop_color, loop_linestyle, sig_marker)
apply_style()

LOOPS = ["cdr1a", "cdr2a", "cdr3a", "cdr1b", "cdr2b", "cdr3b"]
LL = {"cdr1a": "CDR1\u03b1", "cdr2a": "CDR2\u03b1", "cdr3a": "CDR3\u03b1",
      "cdr1b": "CDR1\u03b2", "cdr2b": "CDR2\u03b2", "cdr3b": "CDR3\u03b2"}
REGIONS = ["alpha1 helix", "alpha2 helix", "other G-domain", "beta-sheet floor"]

z = np.load("results/tables/contact_freq_gdomain182_descriptive_p84v2.npz", allow_pickle=True)
F, Y, S, tcr = z["F"], z["Y"], z["S"], z["tcr_position_labels"]
pos = np.array([str(p) for p in z["mhc_positions"]])
num, ins = z["mhc_imgt_number"].astype(int), np.array([str(i) for i in z["mhc_insertion_code"]])
assert bool(z["excludes_alpha3"]) and F.shape == (75, 182)
reg = (pd.read_csv("results/tables/gdomain182_region_classification_p84v2.csv")
         .assign(imgt_position=lambda d: d.imgt_position.astype(str))
         .set_index("imgt_position")["region"].reindex(pos).values)
assert not pd.isna(reg).any()
MIN_ROW_COMPLEXES = 5
nres_row = (S == 1).any(axis=2).sum(axis=0)
rare = nres_row < MIN_ROW_COMPLEXES
assert [str(t) for t in tcr[rare]] == ["cdr3b:111B", "cdr3b:112C", "cdr3b:112B"], tcr[rare]
loop_of = np.array([str(t).split(":")[0] for t in tcr])
runs = []
for i, L in enumerate(loop_of):
    if runs and runs[-1][0] == L: runs[-1][2] = i + 1
    else: runs.append([L, i, i + 1])
assert [r[0] for r in runs] == LOOPS
bnd = [r[1] for r in runs] + [runs[-1][2]]

w82 = ((num >= 50) & (num <= 90)) | ((num >= 1050) & (num <= 1090))
a1, a2 = np.where(w82 & (num < 1000))[0], np.where(w82 & (num >= 1000))[0]
assert (len(a1), len(a2)) == (41, 43)
helix = np.isin(reg, ["alpha1 helix", "alpha2 helix"])
diag = pd.read_csv("results/tables/germline_diagonal_p84v2.csv").set_index("loop")
Fw = F[:, np.r_[a1, a2]]
vmax = float(Fw[~rare].max())                      # hottest pair over rows resolved in >= 5 complexes
vmax_all = float(Fw.max())
hp = np.unravel_index(np.argmax(np.where(rare[:, None], -1, Fw)), Fw.shape)
hp_all = np.unravel_index(np.argmax(Fw), Fw.shape)
wpos = pos[np.r_[a1, a2]]
Fm = np.ma.masked_array(F, mask=np.repeat(rare[:, None], F.shape[1], axis=1))
CM = CMAP_SEQ.copy(); CM.set_bad("white")

fig = plt.figure(figsize=(PAGE_WIDTH_IN, 8.6), constrained_layout=True)
# top room for the panel letters (style.py); this figure leaves constrained layout before the letters are placed
try: fig.get_layout_engine().set(h_pad=0.10, hspace=0.05, rect=(0, 0, 1, 1 - (PANEL_LETTER_OFFSET_IN[1] + PANEL_LETTER_HEIGHT_IN) / 8.6))
except Exception: pass
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.62])

# ---- A -------------------------------------------------------------------
gsa = gs[0, :].subgridspec(2, 3, height_ratios=[0.045, 1], width_ratios=[41, 43, 1.3],
                           hspace=0.02, wspace=0.03)
axA, axStrip = [], []
for k, (cols, ttl, side) in enumerate([(a1, "\u03b11 helix", "alpha1"),
                                       (a2, "\u03b12 helix", "alpha2")]):
    hc = REGION_COLORS[f"{side} helix"]
    axs = fig.add_subplot(gsa[0, k]); axStrip.append(axs)
    axs.imshow(np.array([[to_rgb(hc) if helix[c] else to_rgb(REGION_COLORS["other G-domain"])
                          for c in cols]]), aspect="auto", interpolation="nearest")
    axs.set_xticks([]); axs.set_yticks([])
    for sp in axs.spines.values(): sp.set_visible(False)
    axs.set_title(ttl, fontsize=7.5, pad=3)

    ax = fig.add_subplot(gsa[1, k]); axA.append(ax)
    # sqrt norm: vmax is set by one cell (CDR1b x 76 = 0.667 over rows resolved in >= 5
    # complexes), so a linear scale leaves the rest of the panel near-white. Colorbar ticks still
    # show true frequencies. Rare rows are masked (white) and hatched.
    im = ax.imshow(Fm[:, cols], aspect="auto", cmap=CM,
                   norm=PowerNorm(0.5, vmin=0, vmax=vmax), interpolation="nearest")
    for rr in np.where(rare)[0]:
        ax.add_patch(Rectangle((-0.5, rr - 0.5), len(cols), 1, fill=False, hatch="//////",
                               edgecolor="0.55", lw=0, zorder=4))
    labs = [f"{num[c] % 1000}{ins[c]}" for c in cols]
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([l if j % 2 == 0 else "" for j, l in enumerate(labs)],
                       rotation=90, fontsize=5.4)
    ax.tick_params(axis="x", length=1.5, pad=1)
    for b in bnd[1:-1]: ax.axhline(b - 0.5, color="0.25", lw=0.6)
    if k == 0:
        ax.set_yticks([(bnd[i] + bnd[i + 1]) / 2 - 0.5 for i in range(6)])
        ax.set_yticklabels([LL[L] for L in LOOPS], fontsize=7)
        ax.set_xlabel("IMGT position", fontsize=7.5)
    else:
        ax.set_yticks([])
        ax.set_xlabel("IMGT position (\u03b12 shown without 1000 offset)", fontsize=7.5)
    for i, L in enumerate(LOOPS):
        if diag.loc[L, "prefers"] == side:
            # inset so boxes on adjacent loops stay visually separate
            ax.add_patch(Rectangle((-0.35, bnd[i] - 0.35), len(cols) - 0.3,
                                   (bnd[i + 1] - bnd[i]) - 0.3, fill=False,
                                   edgecolor=loop_color(L), lw=1.2,
                                   ls=loop_linestyle(L), zorder=5, clip_on=False))
cb = fig.colorbar(im, cax=fig.add_subplot(gsa[1, 2]))
cb.set_label("Contact frequency", fontsize=7); cb.ax.tick_params(labelsize=6)

# ---- B / C ---------------------------------------------------------------
def autocrop(path, pad=10):
    im = mpimg.imread(path)
    if im.ndim == 3 and im.shape[2] == 4 and (im[..., 3] < 0.99).any():
        mask = im[..., 3] > 0.02
    else:
        mask = (im[..., :3] < 0.97).any(axis=2)
    ys, xs = np.where(mask)
    return im[max(ys.min() - pad, 0):ys.max() + pad, max(xs.min() - pad, 0):xs.max() + pad]

srcs, axBC = [], []
for k, (kw, ttl) in enumerate([("footprint", "MHC platform by contact frequency (relative)"),
                               ("complex", "TCR\u2013pMHC complex with CDR loops")]):
    ax = fig.add_subplot(gs[1, k]); axBC.append(ax)
    hits = sorted(glob.glob(f"figures/pymol/p84v2_panel*{kw}*.png"))
    if hits:
        ax.imshow(autocrop(hits[0])); srcs.append((kw, hits[0]))
    else:
        ax.text(0.5, 0.5, "PyMOL render pending (55_pymol_fig5_p84v2)", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="0.5"); srcs.append((kw, "MISSING"))
    ax.set_anchor("N"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(ttl, fontsize=7.5, pad=3)
axBC[0].text(0.5, -0.02, f"sticks: > 0.30 of maximum; TCR rows resolved in < {MIN_ROW_COMPLEXES} complexes excluded",
             transform=axBC[0].transAxes, ha="center", va="top", fontsize=6, color="0.3")
axBC[1].legend(handles=[Line2D([], [], marker="s", ls="", ms=5, mec="0.4", mew=0.3,
                               mfc=MOLECULE_COLORS[m], label=lab)
                        for m, lab in (("TCR alpha", "TCR\u03b1"), ("TCR beta", "TCR\u03b2"),
                                       ("peptide", "peptide"), ("MHC", "MHC"))],
               loc="lower center", bbox_to_anchor=(0.5, -0.06), ncol=4, frameon=False,
               fontsize=6, handletextpad=0.3, columnspacing=1.0)

# ---- D -------------------------------------------------------------------
axD = fig.add_subplot(gs[2, 0])
M = np.zeros((6, 4)); K = np.zeros((6, 4), int)
for i, L in enumerate(LOOPS):
    li = loop_of == L
    for j, r in enumerate(REGIONS):
        m = reg == r
        n = int((S[:, li][:, :, m] == 1).sum()); K[i, j] = int(Y[:, li][:, :, m].sum())
        M[i, j] = K[i, j] / n if n else np.nan
im2 = axD.imshow(M * 100, cmap=CMAP_SEQ, aspect="auto", vmin=0)
for i in range(6):
    for j in range(4):
        v = M[i, j] * 100
        axD.text(j, i, "0" if K[i, j] == 0 else f"{v:.2f}", ha="center", va="center",
                 fontsize=6, color="white" if v > np.nanmax(M) * 100 * 0.55 else "0.15")
axD.set_xticks(range(4))
axD.set_xticklabels(["\u03b11\nhelix", "\u03b12\nhelix", "other\nG-domain", "\u03b2-sheet\nfloor"],
                    fontsize=6.5)
axD.set_yticks(range(6)); axD.set_yticklabels([LL[L] for L in LOOPS], fontsize=7)
floor = reg == "beta-sheet floor"
floor_k, floor_n = int(Y[:, :, floor].sum()), int((S[:, :, floor] == 1).sum())
axD.set_xlabel(f"\u03b2-sheet floor: {floor_k} contacts in {floor_n:,} resolved pairs",
               fontsize=6.5)
cb2 = fig.colorbar(im2, ax=axD, fraction=0.05, pad=0.03)
cb2.set_label("Contact rate (%)", fontsize=6.5); cb2.ax.tick_params(labelsize=6)

# ---- E -------------------------------------------------------------------
axE = fig.add_subplot(gs[2, 1])
yy = np.arange(6)[::-1]
for L, yv in zip(LOOPS, yy):
    r = diag.loc[L]; c = loop_color(L)
    sig = not (r.ci_lo <= 0 <= r.ci_hi)
    axE.plot([r.ci_lo * 100, r.ci_hi * 100], [yv, yv], color=c, lw=1.3,
             ls=loop_linestyle(L), solid_capstyle="butt")
    axE.plot(r.mean_diff * 100, yv, "o", ms=4, zorder=3, **sig_marker(c, sig))
    axE.text(1.0, yv, pstr(r.p_bootstrap), transform=axE.get_yaxis_transform(),
             ha="left", va="center", fontsize=6, color="0.3")
axE.axvline(0, color="0.35", lw=0.7)
axE.set_yticks(yy); axE.set_yticklabels([LL[L] for L in LOOPS], fontsize=7)
axE.set_ylim(-0.7, 5.7)
axE.set_xlabel("Per-cluster rate difference, \u03b11 \u2212 \u03b12 (%)", fontsize=7)
axE.text(0.0, 1.03, "\u2190 prefers \u03b12", transform=axE.transAxes, fontsize=6, color="0.35")
axE.text(1.0, 1.03, "prefers \u03b11 \u2192", transform=axE.transAxes, fontsize=6,
         color="0.35", ha="right")

# ---- key, then panel letters --------------------------------------------
fig.canvas.draw()
try: fig.set_layout_engine("none")
except Exception: fig.set_constrained_layout(False)
rend = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
bb = axA[0].get_tightbbox(rend).transformed(inv)
# the alpha2 axis carries the 61A/72A tick labels, so its x-label sits lower: key below both
key_y = min(a.get_tightbbox(rend).transformed(inv).y0 for a in axA) - 0.004
ktxt = fig.text(bb.x0, key_y,
         "box = preferred helix (E)  |  solid = germline, dashed = CDR3  |  color = TCR chain  |  "
         "strip: DSSP helix (dark) / boundary (light)  |  "
         f"hatched: rows resolved in < {MIN_ROW_COMPLEXES} complexes, off scale",
         fontsize=5.8, color="0.3", va="top", ha="left")
fig.canvas.draw()
kb = ktxt.get_window_extent(rend).transformed(inv)
# B and C were laid out under alpha1's (higher) x-label; lower their tops below the key line
for ax in axBC:
    need = ax.get_tightbbox(rend).transformed(inv).y1 - (kb.y0 - 0.024)
    if need > 0:
        p_ = ax.get_position(); ax.set_position([p_.x0, p_.y0, p_.width, p_.height - need])
fig.canvas.draw()
# one fixed offset for every figure (style.panel_letters): A spans its helix strips, heat maps and colour bar; D its colour bar
panel_letters(fig, {"A": axStrip + axA + [cb.ax], "B": [axBC[0]], "C": [axBC[1]], "D": [axD, cb2.ax], "E": [axE]})

png, pdf = save_fig(fig, "figures/fig5_footprint_p84v2_v1")

# ---- validation ----------------------------------------------------------
tot = int(Y[S == 1].sum())
print(f"\npanel B/C sources: {srcs}")
print(f"window columns a1/a2 = {len(a1)}/{len(a2)}; DSSP helix within = "
      f"{int(helix[a1].sum())}/{int(helix[a2].sum())}")
print(f"contacts inside the 82-window: {int(Y[:, :, w82].sum())} / {tot} "
      f"({100 * Y[:, :, w82].sum() / tot:.1f}%)")
print(f"rows resolved in < {MIN_ROW_COMPLEXES} complexes: "
      + ", ".join(f"{t} ({n})" for t, n in zip(tcr[rare], nres_row[rare])))
print(f"hottest pair, rows >= {MIN_ROW_COMPLEXES}: {tcr[hp[0]]} x {wpos[hp[1]]} = {vmax:.3f} (heatmap vmax)")
print(f"hottest pair, all rows: {tcr[hp_all[0]]} x {wpos[hp_all[1]]} = {vmax_all:.3f} (excluded: rare row)")
print(f"beta floor {floor_k} contacts in {floor_n:,} resolved pairs")
mass = F[~rare][:, :].sum(0); mass = mass / mass.max()
print(f"hot positions (> 0.30 relative, rows >= {MIN_ROW_COMPLEXES}): {list(pos[mass > 0.30])}")
lab_a2 = [f"{num[c] % 1000}{ins[c]}" for j, c in enumerate(a2) if j % 2 == 0]
assert "61A" in lab_a2 and "72A" in lab_a2, lab_a2
print("\npanel E markers (fill = CI excludes zero):")
for L in LOOPS:
    r = diag.loc[L]
    print(f"  {LL[L]:6s} {r.prefers:6s} mean {r.mean_diff*100:+7.3f}% "
          f"[{r.ci_lo*100:+.3f},{r.ci_hi*100:+.3f}] {pstr(r.p_bootstrap):11s} "
          f"{'filled' if not (r.ci_lo <= 0 <= r.ci_hi) else 'HOLLOW'} "
          f"{'solid' if L in GERMLINE_LOOPS else 'dashed'}")

from viz.figstats import pdf_font_check
n3, arial, bases = pdf_font_check(pdf)
print(f"\nPDF font check {pdf}: Type3 = {n3}, ArialMT present = {arial}, BaseFonts = {bases}")
assert n3 == 0 and arial
