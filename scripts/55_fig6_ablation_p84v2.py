#!/usr/bin/env python3
"""Figure 6 (p84v2) - sequence dependence, success and failure (tcr70, pairwise task).

Built from 41b_fig6_ablation_p82v1.py; inputs are the p84v2 tables. Under p84v2 three ablations
survive Holm across the 8 (shuffle_tcr, shuffle_mhc, no_germline; under p82v1 shuffle_mhc and
no_germline). The showcase complexes are re-derived by the same rule (largest and smallest
per-complex delta in percomplex_tcr70_paired_p84v2.tsv): 3qfj and 3qeq, no longer 5yxn and 7n1f.

A  ablation forest: 8 ablations vs the full model, delta = full - ablation (positive = the
   ablation lowers macro-AUPRC). Bars are Bonferroni-simultaneous 99.375 % percentile intervals
   (quantiles 0.003125 / 0.996875 of the SAME 20,000 cluster-bootstrap resamples that give the
   95 % CI and p; author decision, revision 1); fill = that simultaneous interval excludes zero.
   Raw p and Holm-adjusted p across the 8 are printed beside each row; the script checks that the
   fill verdict agrees with Holm < 0.05 for every row.
   Recomputed at 20,000 resamples (seed 42) from the per-complex ablation arrays cached by
   scripts/53_ablations_p84v2.py (results/tables/_partial/p84v2/29b_ablations/*.npz); those
   arrays first reproduce every value in ablations_p84v2.csv exactly at 2,000 resamples/seed 0.
   no_cross and dot_head (interval spans zero) are a finding: cross-attention and the bilinear head
   are not demonstrated necessary.
B  per-complex delta (TCR-IMAP - pair-frequency prior, pairwise AUPRC), tcr70, all 174
   complexes ranked; 3qfj (+0.279) and 3qeq (-0.158) marked; pooled delta at 20k.
C  3qfj, observed vs TCR-IMAP-predicted contact sites (largest gain over the prior)
D  3qeq, same (largest loss to the prior). Panels C/D: renders from
   scripts/55_pymol_fig6_p84v2.py (residues selected by IMGT number + insertion code), tcr70 out-of-fold predictions only; predicted set = top-k
   resolved pairs by TCR-IMAP score, k = that complex's number of observed contacts.

Panel A fill = simultaneous interval excludes zero (sig_marker); hollow means only that. Panel B
has no markers. p displayed only through style.pstr(). Colour: model variants = TCR-IMAP
blue; panel B bars by which method is higher (TCR-IMAP blue / prior grey); C/D from
MOLECULE_COLORS (colour = TCR chain).
"""
import sys, glob
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.lines import Line2D
sys.path.insert(0, "src")
from viz.style import (apply_style, save_fig, print_validation_table, pstr, sig_marker, signed, minus,
                       PAGE_WIDTH_IN, PANEL_LABEL_STYLE as PL, MODEL_LABEL, METHOD_COLORS, panel_letters,
                       MOLECULE_COLORS)
from viz.figstats import (pooled_cluster_bootstrap as pcb, holm, ci_excludes_zero,
                          pdf_font_check, N_BOOT, SEED)

apply_style()
T = "results/tables/"
OUT_STEM = "figures/fig6_ablation_p84v2_v1"
CM, CP = METHOD_COLORS[MODEL_LABEL], METHOD_COLORS["Pair-frequency prior"]
TXT = "0.30"
SUCCESS, FAILURE = "3qfj", "3qeq"     # re-derived from p84v2, asserted below

# ------------------------------------------------------------------ ablations
ABL = [("input masked", "no_mhc", "MHC residue identity masked"),
       ("input masked", "no_tcr", "TCR residue identity masked"),
       ("shuffled", "shuffle_tcr", "TCR residues shuffled within loops"),
       ("shuffled", "shuffle_mhc", "MHC residues shuffled across positions"),
       ("loops removed", "no_germline", "CDR1/2 removed"),
       ("loops removed", "no_cdr3", "CDR3 removed"),
       ("architecture", "no_cross", "cross-attention removed"),
       ("architecture", "dot_head", "dot-product head (no bilinear map)")]
rec = pd.read_csv(T + "ablations_p84v2.csv").set_index("ablation")
CACHE = T + "_partial/p84v2/29b_ablations/"
full = np.load(CACHE + "full.npz")
pw = np.load(T + "pairwise_ours_p84v2.npz", allow_pickle=True)
md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
assert (full["idx"] == np.arange(174)).all() and np.allclose(full["auprc"], pw["auprc"], equal_nan=True)
cl = md.loc[pw["pdb"], "clu70"].values
assert (cl == pw["clu70"]).all() and len(np.unique(cl)) == 69
assert abs(round(float(np.nanmean(full["auprc"])), 4) - round(rec.loc["full", "auprc"], 4)) < 1e-9

SIM_Q = 0.05 / 8 / 2          # 0.003125 -> 99.375 % simultaneous interval
res, val = {}, []
for _, key, _ in ABL:
    a = np.load(CACHE + f"{key}.npz")
    assert (a["idx"] == full["idx"]).all()
    r2 = pcb(full["auprc"], a["auprc"], cl, 2000, 0)
    r20 = pcb(full["auprc"], a["auprc"], cl, N_BOOT, SEED, return_draws=True)
    # Bonferroni-simultaneous 99.375 % interval (alpha 0.05 / 8 ablations) from the same draws
    sim_lo, sim_hi = (float(x) for x in np.quantile(r20.pop("draws"), [SIM_Q, 1 - SIM_Q]))
    assert sim_lo <= r20["lo"] and sim_hi >= r20["hi"]
    r20["sim_lo"], r20["sim_hi"] = sim_lo, sim_hi
    rr = rec.loc[key]
    rep = (abs(r2["diff"] - rr["diff"]) < 1e-12 and abs(r2["lo"] - rr.lo) < 1e-12
           and abs(r2["hi"] - rr.hi) < 1e-12 and abs(r2["p"] - rr.p) < 1e-12)
    res[key] = dict(r20=r20, rec=rr, rep=rep, mean=float(np.nanmean(a["auprc"])))
p20 = np.array([res[k]["r20"]["p"] for _, k, _ in ABL])
p_holm = holm(p20)
p_holm_rec = holm(np.array([res[k]["rec"].p for _, k, _ in ABL]))
for (g, k, lab), ph, phr in zip(ABL, p_holm, p_holm_rec):
    res[k]["holm"], res[k]["holm_rec"] = float(ph), float(phr)
    r = res[k]["r20"]
    res[k]["fill"] = ci_excludes_zero(r["sim_lo"], r["sim_hi"])
    res[k]["holm_sig"] = res[k]["holm"] < 0.05
    res[k]["agree"] = res[k]["fill"] == res[k]["holm_sig"]


def ptab(p, n):
    return f"floor ({pstr(p)})" if p <= 1 / n + 1e-12 else f"{p:.4f}"


# ------------------------------------------------------------------ per-complex delta
pr = pd.read_csv(T + "percomplex_tcr70_paired_p84v2.tsv", sep="\t")
assert (pr.pdb.values == pw["pdb"]).all() and np.allclose(pr.grip.values, pw["auprc"])
delta = pr.delta.values
pooled = pcb(pr.grip.values, pr.prior.values, cl, N_BOOT, SEED)
wins = int((delta > 0).sum())
ordr = np.argsort(-delta)
assert pr.pdb.values[ordr[0]] == SUCCESS and pr.pdb.values[ordr[-1]] == FAILURE

# ------------------------------------------------------------------ PyMOL renders + sites
sites = pd.read_csv("figures/pymol/fig6_sites_p84v2.tsv", sep="\t")
SUM = {}
for pdb in (SUCCESS, FAILURE):
    o = sites[(sites.pdb == pdb) & (sites.set == "observed")]
    p = sites[(sites.pdb == pdb) & (sites.set == "predicted")]
    k = len(o)
    assert len(p) == k and int(o.observed.sum()) == k and int(p["rank"].max()) == k
    SUM[pdb] = dict(k=k, hit=int(p.observed.sum()),
                    mhc_o=o.mhc_imgt.nunique(), mhc_p=p.mhc_imgt.nunique(),
                    mhc_shared=len(set(o.mhc_imgt) & set(p.mhc_imgt)),
                    row=pr[pr.pdb == pdb].iloc[0])
    assert SUM[pdb]["k"] == int(SUM[pdb]["row"].npos)


def autocrop(path, pad=12):
    im = mpimg.imread(path)
    mask = im[..., 3] > 0.02 if (im.ndim == 3 and im.shape[2] == 4) else (im[..., :3] < 0.97).any(axis=2)
    ys, xs = np.where(mask)
    im = im[max(ys.min() - pad, 0):ys.max() + pad, max(xs.min() - pad, 0):xs.max() + pad]
    if im.shape[2] == 4:   # composite onto white so the vector PDF has no transparency group
        a = im[..., 3:4]; im = im[..., :3] * a + (1 - a)
    return im


# ------------------------------------------------------------------ figure
FW, FH = PAGE_WIDTH_IN, 6.45
fig = plt.figure(figsize=(FW, FH))
def axin(l, b, w, h): return fig.add_axes([l / FW, b / FH, w / FW, h / FH])
PANELS = {}   # {letter: panel content}; style.panel_letters places every letter at one fixed offset

# --- A -----------------------------------------------------------------------------------
ax = axin(1.62, 3.45, 1.50, 2.55)
ys = []
y = 0.0; prev = None
for g, k, lab in ABL:
    if prev is not None and g != prev: y += 0.55
    ys.append(y); prev = g; y += 1
ys = -np.array(ys)
for (g, k, lab), yv in zip(ABL, ys):
    r = res[k]["r20"]
    ax.plot([r["sim_lo"], r["sim_hi"]], [yv, yv], color=CM, lw=1.1, solid_capstyle="butt")
    ax.plot(r["diff"], yv, "o", ms=4.6, zorder=3, **sig_marker(CM, res[k]["fill"]))
    ax.text(1.03, yv, pstr(r["p"]), transform=ax.get_yaxis_transform(), ha="left", va="center",
            fontsize=5.8, color=TXT)
    ax.text(1.36, yv, pstr(res[k]["holm"]), transform=ax.get_yaxis_transform(), ha="left",
            va="center", fontsize=5.8, color=TXT)
ax.text(1.03, 1.01, "raw", transform=ax.transAxes, ha="left", va="bottom", fontsize=6, color="0.2")
ax.text(1.36, 1.01, "Holm (8)", transform=ax.transAxes, ha="left", va="bottom", fontsize=6, color="0.2")
ax.axvline(0, ls="--", lw=0.8, color="0.4")
ax.set_yticks(ys); ax.set_yticklabels([lab for _, _, lab in ABL], fontsize=6.2)
ax.set_ylim(ys.min() - 0.7, 0.7)
xmin = min(res[k]["r20"]["sim_lo"] for _, k, _ in ABL)
xmax = max(res[k]["r20"]["sim_hi"] for _, k, _ in ABL)
ax.set_xlim(np.floor((xmin - 0.004) * 100) / 100, np.ceil((xmax + 0.004) * 100) / 100)
ax.set_xlabel("Full − ablated (pairwise macro-AUPRC)")
ax.text(0.5, -0.165, "bars: Bonferroni-simultaneous 99.375% intervals (8 ablations);\n"
        "filled = interval excludes 0", transform=ax.transAxes, ha="center", va="top",
        fontsize=6, color=TXT, linespacing=1.25)
yarch = ys[-2:].mean()
xarch = max(res[k]["r20"]["sim_hi"] for k in ("no_cross", "dot_head"))
ax.annotate("interval spans 0:\nnot shown\nnecessary", xy=(xarch + 0.002, yarch),
            xytext=(xarch + 0.010, yarch),
            fontsize=5.6, color=CM, va="center", ha="left", linespacing=1.2,
            arrowprops=dict(arrowstyle="-[, widthB=1.6, lengthB=0.3", lw=0.7, color=CM))
ax.set_title(minus(f"Ablations, tcr70 (full {MODEL_LABEL} {rec.loc['full', 'auprc']:.4f})"), fontsize=8, pad=11)
PANELS["A"] = [ax]

# --- B -----------------------------------------------------------------------------------
ax = axin(4.88, 3.45, 2.17, 2.55)
d_sorted = delta[ordr]
cols = [CM if v > 0 else CP for v in d_sorted]
ax.bar(np.arange(len(d_sorted)), d_sorted, width=1.0, color=cols, linewidth=0)
ax.axhline(0, color="0.2", lw=0.6)
ax.axhline(pooled["diff"], ls="--", lw=0.8, color=CM)
for pdb, ha, dy in ((SUCCESS, "left", -0.02), (FAILURE, "right", 0.02)):
    j = int(np.where(pr.pdb.values[ordr] == pdb)[0][0]); v = d_sorted[j]
    ax.annotate(f"{pdb}  {signed(v, 4)}", xy=(j, v), xytext=(j + (14 if ha == "left" else -14), v + dy),
                fontsize=6, ha=ha, va="center", color="0.15",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.3", shrinkA=0, shrinkB=1))
ax.text(0.98, 0.97, minus(f"{MODEL_LABEL} higher on {wins}/174 ({100 * wins / 174:.1f}%)\n")+
        f"pooled Δ = {signed(pooled['diff'], 4)} (dashed)\n95% CI [{signed(pooled['lo'], 4)}, {signed(pooled['hi'], 4)}]\n"
        f"{pstr(pooled['p'])}", transform=ax.transAxes, ha="right", va="top", fontsize=6,
        color=TXT, linespacing=1.35)
ax.set_xlim(-2, len(d_sorted) + 1); ax.set_xticks([0, 50, 100, 150])
ax.set_xlabel("Complexes, ranked by Δ")
ax.set_ylabel(f"Δ pairwise AUPRC\n({MODEL_LABEL} − pair-frequency prior)", linespacing=1.2)
ax.set_ylim(-0.30, 0.30)
ax.set_title("Per-complex difference, tcr70", fontsize=8, pad=11)
PANELS["B"] = [ax]

# --- C / D -------------------------------------------------------------------------------
IMW, IMH, GAP = 1.58, 1.62, 0.08
x0 = 0.18
xs = [x0, x0 + IMW + GAP, x0 + 2 * (IMW + GAP) + 0.20, x0 + 3 * (IMW + GAP) + 0.20]
srcs, imgax = [], []
for n, (pdb, kind) in enumerate([(SUCCESS, "observed"), (SUCCESS, "predicted"),
                                 (FAILURE, "observed"), (FAILURE, "predicted")]):
    ax = axin(xs[n], 0.55, IMW, IMH); imgax.append(ax)
    path = f"figures/pymol/p84v2_fig6_{pdb}_{kind}.png"
    ax.imshow(autocrop(path)); srcs.append(path)
    ax.set_anchor("N"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    s = SUM[pdb]
    ax.set_title("observed contacts (k = %d)" % s["k"] if kind == "observed"
                 else f"{MODEL_LABEL} top {s['k']} pairs", fontsize=6.8, pad=2)
    if n == 0:   # orientation, identical for every render (32c derives the frame per structure)
        ax.text(0.02, 0.98, "α1", transform=ax.transAxes, fontsize=6.5, color="0.25", va="top")
        ax.text(0.02, 0.02, "α2", transform=ax.transAxes, fontsize=6.5, color="0.25", va="bottom")
        ax.annotate("", xy=(0.98, 0.04), xytext=(0.74, 0.04), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", lw=0.7, color="0.35"))
        ax.text(0.86, 0.075, "peptide N→C", transform=ax.transAxes, fontsize=5.2,
                color="0.35", ha="center", va="bottom")
for (pdb, xl, lab, tag), pair in zip(((SUCCESS, xs[0], "C", "largest gain over the prior"),
                                     (FAILURE, xs[2], "D", "largest loss to the prior")), (imgax[:2], imgax[2:])):
    r = SUM[pdb]["row"]
    s = SUM[pdb]
    hdr = fig.text((xl + IMW + GAP / 2) / FW, 2.42 / FH,
             minus(f"{pdb}: {tag}\n{MODEL_LABEL} {r.grip:.3f} vs prior {r.prior:.3f}") + f" (Δ = {signed(r.delta, 4)})",
             ha="center", va="bottom", fontsize=6.6, color="0.15", linespacing=1.3)
    ftr = fig.text((xl + IMW + GAP / 2) / FW, 0.50 / FH,
             minus(f"{s['hit']} of {s['k']} predicted pairs observed; MHC sites shared {s['mhc_shared']} of {s['mhc_o']}"),
             ha="center", va="top", fontsize=5.8, color=TXT)
    PANELS[lab] = list(pair) + [hdr, ftr]
handles = [Line2D([], [], marker="s", ls="", ms=5.5, mfc=MOLECULE_COLORS["TCR alpha"], mec="none",
                  label="contacted by TCRα"),
           Line2D([], [], marker="s", ls="", ms=5.5, mfc=MOLECULE_COLORS["TCR beta"], mec="none",
                  label="contacted by TCRβ"),
           Line2D([], [], marker="o", ls="", ms=4.5, mfc="0.45", mec="none",
                  label="CDR contact residue (Cα, chain color)"),
           Line2D([], [], marker="s", ls="", ms=5.5, mfc=MOLECULE_COLORS["MHC"], mec="none",
                  label="MHC α1/α2 surface"),
           Line2D([], [], marker="s", ls="", ms=5.5, mfc=MOLECULE_COLORS["peptide"], mec="none",
                  label="peptide")]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=5, frameon=False,
           fontsize=5.8, handletextpad=0.3, columnspacing=1.1)

panel_letters(fig, PANELS)          # one fixed offset for every figure (style.py)
png, pdf = save_fig(fig, OUT_STEM)
plt.close(fig)

# ------------------------------------------------------------------ validation
rows = []
for g, k, lab in ABL:
    r, rr = res[k]["r20"], res[k]["rec"]
    rows.append(dict(abl=k, mean=f"{res[k]['mean']:.4f}/{rr.auprc:.4f}",
                     rec=f"{rr['diff']:+.4f} [{rr.lo:+.4f},{rr.hi:+.4f}] {ptab(rr.p, 2000)}",
                     rep="exact" if res[k]["rep"] else "MISMATCH",
                     r20=f"{r['diff']:+.4f} [{r['lo']:+.4f},{r['hi']:+.4f}] {ptab(r['p'], N_BOOT)}",
                     holm20=pstr(res[k]["holm"]), holmrec=pstr(res[k]["holm_rec"]),
                     sim=f"[{r['sim_lo']:+.4f},{r['sim_hi']:+.4f}]",
                     shown=f"{pstr(r['p'])} | {pstr(res[k]['holm'])}",
                     ci95="excl 0" if ci_excludes_zero(r["lo"], r["hi"]) else "incl 0",
                     fill="filled" if res[k]["fill"] else "hollow",
                     holmv="< 0.05" if res[k]["holm_sig"] else ">= 0.05",
                     agree="agree" if res[k]["agree"] else "DISAGREE"))
print("\n--- panel A: ablations (recorded = ablations_p84v2.csv, 2k; recomputed 20k seed 42) ---")
print_validation_table(rows, [("abl", "ablation", None), ("mean", "AUPRC calc/rec", None),
                              ("rec", "recorded 2k", None), ("rep", "2k repro", None),
                              ("r20", "recomputed 20k", None), ("holmrec", "Holm on rec", None),
                              ("holm20", "Holm on 20k", None), ("shown", "displayed raw | Holm", None),
                              ("ci95", "95% CI", None), ("sim", "99.375% simultaneous (20k)", None),
                              ("fill", "marker (sim.)", None), ("holmv", "Holm verdict", None),
                              ("agree", "fill vs Holm", None)])
dis = [k for _, k, _ in ABL if not res[k]["agree"]]
print(f"fill (simultaneous interval excludes 0) vs Holm < 0.05: "
      f"{'all 8 rows agree' if not dis else 'DISAGREE in ' + ', '.join(dis)}")
print("Holm-adjusted p (20k), unrounded: " + ", ".join(f"{k} {res[k]['holm']:.4f}" for _, k, _ in ABL))
print(f"\n--- panel B: percomplex_tcr70_paired_p84v2.tsv ---\n"
      f"pooled delta {pooled['diff']:+.4f} [{pooled['lo']:+.4f},{pooled['hi']:+.4f}] {ptab(pooled['p'], N_BOOT)} "
      f"(shown {pstr(pooled['p'])}); wins {wins}/174; max {pr.pdb.values[ordr[0]]} {d_sorted[0]:+.4f}; "
      f"min {pr.pdb.values[ordr[-1]]} {d_sorted[-1]:+.4f}")
print("\n--- panels C/D: figures/pymol/fig6_sites_p84v2.tsv (55_pymol_fig6_p84v2) ---")
for pdb in (SUCCESS, FAILURE):
    s = SUM[pdb]; r = s["row"]
    print(f"  {pdb}: k = {s['k']} observed contacts = npos {int(r.npos)}; top-{s['k']} predicted pairs observed "
          f"{s['hit']}/{s['k']}; MHC sites obs {s['mhc_o']} pred {s['mhc_p']} shared {s['mhc_shared']}; "
          f"TCR-IMAP {r.grip:.4f} prior {r.prior:.4f} delta {r.delta:+.4f}")
print(f"  renders: {srcs}")

# cross-check against the step-1 table (scripts/55_stats_p84v2.py, 20,000 resamples, seed 42)
S20 = pd.read_csv("results/review/figures_p84v2/stats_20k_seed42.tsv", sep="\t").set_index(["group", "comparison"])
chk = []
for _, k, _ in ABL:
    s, r = S20.loc[("ablation (full - ablated)", k)], res[k]["r20"]
    chk.append((f"ablation {k}: diff/CI/p/Holm/simultaneous",
                all(abs(r[a] - s[b]) < 1e-12 for a, b in (("diff", "diff"), ("lo", "lo"), ("hi", "hi"), ("p", "p"),
                                                         ("sim_lo", "sim_lo"), ("sim_hi", "sim_hi")))
                and abs(res[k]["holm"] - s.p_holm) < 1e-12 and res[k]["fill"] == bool(s.sim_excludes_zero)))
s = S20.loc[("pairwise", "tcr70: TCR-IMAP vs pair-frequency prior")]
chk.append(("pairwise pooled vs prior", all(abs(pooled[a] - s[a]) < 1e-12 for a in ("diff", "lo", "hi", "p"))
            and f"{wins}/174" == s.wins))
surv = [k for _, k, _ in ABL if res[k]["holm"] < 0.05]
chk.append(("Holm survivors == shuffle_tcr, shuffle_mhc, no_germline", surv == ["shuffle_tcr", "shuffle_mhc", "no_germline"]))
chk.append(("all 8 fills agree with Holm", not dis))
chk.append(("extremes re-derived: max 3qfj, min 3qeq", pr.pdb.values[ordr[0]] == "3qfj" and pr.pdb.values[ordr[-1]] == "3qeq"))
print("\nstep-1 (stats_20k_seed42.tsv) cross-check:")
for name, ok in chk:
    print(f"  {'OK      ' if ok else 'MISMATCH'} {name}")
assert all(ok for _, ok in chk)
for old in ("5yxn", "7n1f"):
    j = int(np.where(pr.pdb.values[ordr] == old)[0][0])
    print(f"  p82v1 showcase {old}: p84v2 rank {j + 1}/174, delta {d_sorted[j]:+.4f}")
print("  top 3 / bottom 3: " + ", ".join(f"{pr.pdb.values[ordr[j]]} {d_sorted[j]:+.4f} (clu70 {pr.clu70.values[ordr[j]]})"
                                         for j in (0, 1, 2, -3, -2, -1)))
n3, arial, bases = pdf_font_check(pdf)
print(f"\nPDF font check {pdf}: Type3 = {n3}, ArialMT present = {arial}, BaseFonts = {bases}")
assert n3 == 0 and arial
