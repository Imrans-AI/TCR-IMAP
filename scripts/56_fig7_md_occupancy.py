#!/usr/bin/env python3
"""Figure 7 - MD occupancy of crystal-negative pairs against TCR-IMAP and the pair-frequency prior.

Reads the tables written by the single authorised run (results/tables/md_*.csv; manifest commit 4d53879).
Nothing here recomputes any occupancy, rho or p.

A  Spearman rho by crystal-frame distance stratum, pooled over the four systems, and over all pairs; TCR-IMAP
   and prior side by side; per-system values overlaid. All-pairs group set apart so it is not read as a stratum.
   Fill = approximate within-system permutation p < 0.05 (the established fill-means-significance rule).
B  Score vs system occupancy (mean of 3 replicates) for every crystal-negative pair, pooled over the four
   systems, coloured by stratum; one small multiple per score. Non-zero occupancies on a log axis, the zero
   occupancies in a strip beneath, so the sparsity is visible.
C  Top-k: fraction of the top-k crystal-negative pairs with occupancy >= 0.25 against the base rate, per system,
   both scores; k = crystal-positive count (predeclared). Ties at k marked.
D  Check 1: crystal-positive occupancy per replicate and system, with the docked-frame fractions under both
   readings annotated.
E  Backbone RMSD (protein Ca vs the energy-minimised structure), 20-200 ns, per system and replicate, from the
   characterisation tables recomputed with periodic-boundary handling (md_characterisation_manifest.csv); dashed line
   and label = system mean (mean of the three replicate means, equal frame counts); docked-frame fraction (count
   reading) per replicate, from panel D's table, so the TCR's partial disengagement in 5yxn can be read against it.
   Evidence of simulation stability; descriptive.

Panel B needs two things no md_* table holds: each pair's score, and its distance stratum. Scores are read from
results/tables/pairmaps_tcr70_p84v2.npz after its sha256 is checked against the run manifest's scores_sha256, i.e.
the exact file the run used. Strata use the analysis's own crystal_distance_matrix and stats.stratum_masks. System
occupancy per pair is the predeclared mean of the three replicate values in md_occupancy_pairs.csv. Before anything
is drawn, the resulting per-pair arrays must reproduce the tables' n, n_distinct, n_nonzero and n_ge_025 for every
system x stratum, every pooled stratum, every system and the pooled set (20 groups); otherwise the script stops.
No rho is recomputed, even as a check.

Encoding. TCR-IMAP blue, pair-frequency prior dark grey (METHOD_COLORS). Strata: STRATUM_COLORS (brown ramp, dark =
near). Fill carries significance only in A; C and D are descriptive and drawn with lines, ticks and crosses, which
carry no fill.
"""
import hashlib
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, "."); sys.path.insert(0, "src"); sys.path.insert(0, "md/analysis")
from viz.style import (apply_style, save_fig, print_validation_table, PAGE_WIDTH_IN, PANEL_LABEL_STYLE as PL,
                       MODEL_LABEL, METHOD_COLORS, STRATUM_COLORS, REPLICATE_COLORS, panel_letters, REFERENCE_LINE_STYLE)
from viz.figstats import pdf_font_check
import stats as st                                     # the analysis's own functions (md/analysis/stats.py)
from mdmap_p84v2 import crystal_distance_matrix

apply_style()
T = "results/tables/"
OUT_STEM = "figures/fig7_md_occupancy_v1"
SYSTEMS = ["5yxn", "7n1f", "3vxs", "8eo8"]             # queue order: high, low, high, low Delta (p82v1 selection)
SCORE_KEYS = [("TCR-IMAP", "pred_model"), ("pair-frequency prior", "pred_freq")]
SLAB = {"TCR-IMAP": MODEL_LABEL, "pair-frequency prior": "Pair-frequency prior"}
COL = {"TCR-IMAP": METHOD_COLORS[MODEL_LABEL], "pair-frequency prior": METHOD_COLORS["Pair-frequency prior"]}
STRATA = ["4.5-6 A", "6-8 A", "> 8 A"]
SDISP = {"4.5-6 A": "4.5–6 Å", "6-8 A": "6–8 Å", "> 8 A": "> 8 Å"}
TXT = "0.30"
VAL = []


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def val(panel, item, shown, source, cell):
    VAL.append(dict(panel=panel, item=item, shown=shown, source=source, cell=cell))


# ------------------------------------------------------------------ the run's tables, and provenance
# float_precision="round_trip": pandas' default fast parser misreads ~1 in 70 occupancy values by one ulp, which
# breaks ties that existed in the run's memory; round-trip parsing returns exactly k/9001 for every value
RD = dict(float_precision="round_trip")
man = pd.read_csv(T + "md_occupancy_manifest.csv").set_index("item")["value"]
assert man["md_analysis_tree_clean"] == "True"
assert man["md_occupancy_pairs_sha256_before_scores"] == sha(T + "md_occupancy_pairs.csv"), "pairs table changed"
# The characterisation tables were recomputed with periodic-boundary handling after the run (criteria document,
# 24 October 2026); the run manifest, left untouched, describes the superseded summary, which is kept as a backup.
# The current tables are checked against md_characterisation_manifest.csv in panel E.
assert man["md_characterisation_summary_sha256_before_scores"] == sha(
    "results/review/characterisation_fix/backup_superseded/md_characterisation_summary.csv")
rho = pd.read_csv(T + "md_occupancy_rho.csv", **RD)
strata = pd.read_csv(T + "md_occupancy_strata.csv", **RD)
topk = pd.read_csv(T + "md_occupancy_topk.csv", **RD)
chk1 = pd.read_csv(T + "md_occupancy_check1_positive.csv", dtype={"rep": str}, **RD)
dock = pd.read_csv(T + "md_occupancy_docked_frames.csv", **RD)
pairs = pd.read_csv(T + "md_occupancy_pairs.csv", dtype={"tcr_slot": str, "mhc_position": str}, **RD)
print(f"run commit {man['git_commit'][:7]}; md/analysis clean {man['md_analysis_tree_clean']}; "
      f"phase-A and A2 hashes match the files on disk")

SCORES = T + "pairmaps_tcr70_p84v2.npz"
LABELS = "data/processed/pair_maps_p84v2.npz"
assert sha(SCORES) == man["scores_sha256"], "score file is not the one the run used"
assert sha(LABELS) == man["labels_sha256"], "label file is not the one the run used"
pm = np.load(SCORES, allow_pickle=True); lab = np.load(LABELS, allow_pickle=True)
lpdb = [str(x) for x in lab["pdb"]]
assert [str(x) for x in pm["pdb"]] == lpdb
rowi = {str(x): i for i, x in enumerate(lab["row_label"])}; coli = {str(x): j for j, x in enumerate(lab["col_label"])}
print("score and label files: sha256 equal to the manifest's (the exact inputs of the run)")

# ------------------------------------------------------------------ panel B arrays, and their cross-check
REP = (1, 2, 3)
per = {}
for p in SYSTEMS:
    i = lpdb.index(p)
    Y, S = lab["Y"][i].astype(bool), lab["S"][i].astype(bool)
    assert (pm["Y"][i].astype(bool) == Y).all() and (pm["S"][i].astype(bool) == S).all()
    O = np.full(Y.shape, np.nan)
    acc = np.zeros(Y.shape); seen = np.zeros(Y.shape, int)
    for r in REP:
        sub = pairs[(pairs.pdb == p) & (pairs.rep == r)]
        ti = sub.tcr_slot.map(rowi).values; mj = sub.mhc_position.map(coli).values
        assert not np.isnan(ti.astype(float)).any() and not np.isnan(mj.astype(float)).any()
        assert (sub.crystal_contact.values == Y[ti, mj]).all(), f"{p} rep{r}: crystal_contact differs from Y"
        acc[ti, mj] += sub.occupancy.values; seen[ti, mj] += 1
    assert ((seen == 3) == S).all(), f"{p}: pairs table does not cover exactly the resolved pairs x 3 replicates"
    O[S] = acc[S] / 3.0                                   # predeclared system occupancy = mean of 3 replicates
    d = crystal_distance_matrix(p)
    per[p] = dict(Y=Y, S=S, O=O, D=d, masks=st.stratum_masks(Y, S, d), neg=st.neg_mask(Y, S),
                  sc={name: pm[key][i].astype(float) for name, key in SCORE_KEYS})

KEYS = ["n_distinct_occupancy", "n_nonzero_occupancy", "n_occupancy_ge_025"]
xc_rows, n_bad = [], 0


def xcheck(group, vals, row, ncol):
    global n_bad
    got = dict(n=len(vals), **st.occupancy_counts(vals))
    want = dict(n=int(row[ncol]), **{k: int(row[k]) for k in KEYS})
    ok = got == want
    n_bad += not ok
    xc_rows.append(dict(group=group, n=f"{got['n']}/{want['n']}",
                        distinct=f"{got['n_distinct_occupancy']}/{want['n_distinct_occupancy']}",
                        nonzero=f"{got['n_nonzero_occupancy']}/{want['n_nonzero_occupancy']}",
                        ge025=f"{got['n_occupancy_ge_025']}/{want['n_occupancy_ge_025']}", match="yes" if ok else "NO"))


S1 = strata[strata.score == "TCR-IMAP"]
R1 = rho[(rho.score == "TCR-IMAP")]
for p in SYSTEMS:
    for s_ in STRATA:
        xcheck(f"{p} {s_}", per[p]["O"][per[p]["masks"][s_]],
               S1[(S1.pdb == p) & (S1.stratum == s_)].iloc[0], "n")
for s_ in STRATA:
    xcheck(f"pooled {s_}", np.concatenate([per[p]["O"][per[p]["masks"][s_]] for p in SYSTEMS]),
           S1[(S1.level == "pooled") & (S1.stratum == s_)].iloc[0], "n")
for p in SYSTEMS:
    xcheck(f"{p} all", per[p]["O"][per[p]["neg"]], R1[(R1.level == "system") & (R1.pdb == p)].iloc[0], "n_crystal_negative")
xcheck("pooled all", np.concatenate([per[p]["O"][per[p]["neg"]] for p in SYSTEMS]),
       R1[R1.level == "pooled"].iloc[0], "n_crystal_negative")
print("\nPanel B cross-check - per-pair arrays (this figure) / table cell (run):")
print_validation_table(xc_rows, [("group", "group", None), ("n", "n", None), ("distinct", "n_distinct", None),
                                 ("nonzero", "n_nonzero", None), ("ge025", "n >= 0.25", None), ("match", "match", None)])
assert n_bad == 0, f"{n_bad} groups do not reproduce the tables - STOP"
print(f"all {len(xc_rows)} groups reproduce the tables exactly")

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(PAGE_WIDTH_IN, 7.7))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.15], height_ratios=[1.0, 0.95], hspace=0.78, wspace=0.30,
                      left=0.075, right=0.985, top=0.955, bottom=0.335)
gsE = fig.add_gridspec(1, 4, wspace=0.10, left=0.075, right=0.985, top=0.222, bottom=0.055)

# ---------------- A: rho by stratum, and all pairs
axA = fig.add_subplot(gs[0, 0])
GROUPS = STRATA + ["all"]
XG = {"4.5-6 A": 0.0, "6-8 A": 1.0, "> 8 A": 2.0, "all": 3.55}
W = 0.36
any_ni = False
for g in GROUPS:
    for k, (name, _) in enumerate(SCORE_KEYS):
        x = XG[g] + (k - 0.5) * W
        if g == "all":
            row = rho[(rho.score == name) & (rho.level == "pooled")].iloc[0]
            sysrows = rho[(rho.score == name) & (rho.level == "system")]
            src = "md_occupancy_rho.csv"
        else:
            row = strata[(strata.score == name) & (strata.level == "pooled") & (strata.stratum == g)].iloc[0]
            sysrows = strata[(strata.score == name) & (strata.level == "system") & (strata.stratum == g)]
            src = "md_occupancy_strata.csv"
        sig = bool(row.perm_p_approximate < 0.05)
        # Bars are always solid: the fill rule applies to the single-system markers only (legend). Every pooled
        # p is < 0.05, so this draws exactly what the earlier fill-by-p bars drew; stop if that ever changes.
        assert sig, (name, g, row.perm_p_approximate)
        axA.bar(x, row.rho, W * 0.92, color=COL[name], edgecolor=COL[name], lw=1.0, zorder=2)
        ytop = max(row.rho, sysrows.rho.max())             # above the bar AND its per-system dots
        axA.text(x, ytop + 0.012, f"{row.rho:.2f}", ha="center", va="bottom", fontsize=5.8, color=COL[name], zorder=4)
        val("A", f"{SLAB[name]} {g} pooled rho", f"{row.rho:.2f} (solid bar; p {row.perm_p_approximate:.4f})",
            src, f"score={name}, level=pooled{'' if g == 'all' else ', stratum=' + g}")
        if not bool(row.interpretable):
            any_ni = True
            axA.text(x, row.rho + 0.06, "n.i.", ha="center", fontsize=5.5, color="0.1")
        jit = np.linspace(-W * 0.28, W * 0.28, len(SYSTEMS))
        for jj, p in enumerate(SYSTEMS):
            sr = sysrows[sysrows.pdb == p].iloc[0]
            ssig = bool(sr.perm_p_approximate < 0.05)
            axA.plot(x + jit[jj], sr.rho, "o", ms=2.6, zorder=3, **dict(mec="0.15", mew=0.5,
                     mfc=COL[name] if ssig else "white"))
            if not bool(sr.interpretable):
                any_ni = True
                axA.text(x + jit[jj], sr.rho + 0.03, "n.i.", ha="center", fontsize=5, color="0.1")
            val("A", f"{SLAB[name]} {g} {p}", f"{sr.rho:.3f} ({'filled' if ssig else 'hollow'})", src,
                f"score={name}, level=system, pdb={p}{'' if g == 'all' else ', stratum=' + g}")
    # n and non-zero beneath each group
    if g == "all":
        rr = rho[(rho.score == "TCR-IMAP") & (rho.level == "pooled")].iloc[0]
        n_, nz = int(rr.n_crystal_negative), int(rr.n_nonzero_occupancy)
    else:
        rr = strata[(strata.score == "TCR-IMAP") & (strata.level == "pooled") & (strata.stratum == g)].iloc[0]
        n_, nz = int(rr.n), int(rr.n_nonzero_occupancy)
    axA.text(XG[g], -0.115, f"n = {n_:,}\n{nz} non-zero", ha="center", va="top", fontsize=5.6, color=TXT,
             transform=axA.get_xaxis_transform())
    val("A", f"{g} n / non-zero", f"{n_:,} / {nz}", "md_occupancy_strata.csv" if g != "all" else "md_occupancy_rho.csv",
        "pooled row, n (n_crystal_negative) and n_nonzero_occupancy")
axA.axvline(2.78, color="0.55", lw=0.6, ls=(0, (2, 2)))
axA.set_xticks([XG[g] for g in GROUPS])
axA.set_xticklabels([SDISP[g] for g in STRATA] + ["all pairs"])
axA.text(1.0, 1.0, "crystal-frame distance stratum", transform=axA.get_xaxis_transform(), ha="center",
         va="bottom", fontsize=6.2, color=TXT)
axA.text(XG["all"], 1.0, "pooled", transform=axA.get_xaxis_transform(), ha="center", va="bottom", fontsize=6.2,
         color=TXT)
axA.set_ylabel("Spearman ρ, score vs occupancy\n(crystal-negative pairs)")
axA.set_ylim(0, 0.55); axA.set_xlim(-0.55, 4.05)
axA.legend(handles=[Patch(fc=COL[n], ec=COL[n], label=SLAB[n]) for n, _ in SCORE_KEYS] +
           [Line2D([], [], marker="o", ls="", ms=2.6, mfc="0.5", mec="0.15", mew=0.5,
                   label="single system (filled: approx. p < 0.05)")],     # fill rule = markers only; bars solid
           loc="upper center", bbox_to_anchor=(0.5, -0.27), ncol=3, frameon=False, handlelength=1.0,
           columnspacing=0.9, handletextpad=0.3, fontsize=5.8)

# ---------------- B: score vs occupancy, small multiples
axB_all = []
gsB = gs[0, 1].subgridspec(2, 2, height_ratios=[1.0, 0.24], hspace=0.06, wspace=0.12)
OMIN = 1.0 / (3 * 9001)                                  # smallest non-zero system occupancy possible
rng = np.random.default_rng(0)                           # jitter for the zero strip only (display)
cnt_txt = {}
for c, (name, _) in enumerate(SCORE_KEYS):
    ax = fig.add_subplot(gsB[0, c]); axz = fig.add_subplot(gsB[1, c], sharex=ax); axB_all += [ax, axz]
    for s_ in STRATA[::-1]:                              # far first, near on top
        sc = np.concatenate([per[p]["sc"][name][per[p]["masks"][s_]] for p in SYSTEMS])
        oc = np.concatenate([per[p]["O"][per[p]["masks"][s_]] for p in SYSTEMS])
        nz = oc > 0
        ax.scatter(sc[nz], oc[nz], s=4 if s_ != "> 8 A" else 3.5, c=STRATUM_COLORS[s_], lw=0.25,
                   edgecolors="0.25" if s_ == "> 8 A" else "none", zorder=3 if s_ == "4.5-6 A" else 2)
        axz.scatter(sc[~nz], rng.uniform(-1, 1, (~nz).sum()), s=1.2, c=STRATUM_COLORS[s_], lw=0,
                    alpha=0.6, rasterized=True)
        cnt_txt[s_] = (int(nz.sum()), len(oc))
    ax.set_yscale("log"); ax.set_ylim(OMIN * 0.7, 1.15)
    ax.axhline(0.25, color="0.55", lw=0.6, ls=(0, (3, 2)))
    if c == 0:
        ax.set_ylabel("MD occupancy\n(non-zero, log)")
        axz.set_yticks([0]); axz.set_yticklabels(["0"])
    else:
        ax.tick_params(labelleft=False); axz.set_yticks([])
    axz.set_ylim(-1.6, 1.6)
    ax.tick_params(labelbottom=False)
    axz.spines["left"].set_visible(c == 0)
    ax.set_title(SLAB[name], fontsize=7, color=COL[name], pad=2)
    axz.set_xlabel(f"{'TCR-IMAP score' if name == 'TCR-IMAP' else 'Pair-frequency prior score'}", fontsize=7)
    if c == 0:
        axB0 = ax
for s_ in STRATA:
    val("B", f"{s_} non-zero / pairs (pooled)", f"{cnt_txt[s_][0]} / {cnt_txt[s_][1]:,}", "md_occupancy_strata.csv",
        f"pooled, stratum={s_}: n_nonzero_occupancy / n (cross-checked above)")
_pz = axz.get_position(); _p0 = axB0.get_position()
fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=3, mfc=STRATUM_COLORS[s_],
                           mec="0.25" if s_ == "> 8 A" else STRATUM_COLORS[s_], mew=0.3,
                           label=f"{SDISP[s_]}: {cnt_txt[s_][0]}/{cnt_txt[s_][1]:,} non-zero") for s_ in STRATA] +
                   [Line2D([], [], color="0.55", lw=0.6, ls=(0, (3, 2)), label="occupancy 0.25")],
           loc="upper center", bbox_to_anchor=((_p0.x0 + _pz.x1) / 2, _pz.y0 - 0.058), ncol=2, frameon=False,
           fontsize=5.6, handletextpad=0.3, columnspacing=1.2, handlelength=1.4)

# ---------------- C: top-k against base rate
axC = fig.add_subplot(gs[1, 0])
for i, p in enumerate(SYSTEMS):
    base = topk[(topk.pdb == p)].base_rate.unique()
    assert len(base) == 1
    base = float(base[0])
    axC.plot([i - 0.38, i + 0.38], [base, base], solid_capstyle="butt", zorder=1, **REFERENCE_LINE_STYLE)   # not a method colour
    val("C", f"{p} base rate", f"{base:.5f}", "md_occupancy_topk.csv", f"pdb={p}, base_rate (same for both scores)")
    for k, (name, _) in enumerate(SCORE_KEYS):
        r = topk[(topk.pdb == p) & (topk.score == name)].iloc[0]
        x = i + (k - 0.5) * 0.36
        axC.plot([x, x], [base, r.frac_top_k], color=COL[name], lw=1.4, zorder=2)
        axC.plot([x - 0.1, x + 0.1], [r.frac_top_k, r.frac_top_k], color=COL[name], lw=1.6, solid_capstyle="butt", zorder=3)
        lab_ = f"{int(r.n_top_k_ge)}/{int(r.k)}" + ("†" if bool(r.tie_at_k) else "")
        axC.text(x, r.frac_top_k * 1.25, lab_, ha="center", va="bottom", fontsize=5.6, color=COL[name])
        val("C", f"{p} {SLAB[name]}", f"{r.frac_top_k:.4f} ({lab_})", "md_occupancy_topk.csv",
            f"pdb={p}, score={name}: frac_top_k, n_top_k_ge, k, tie_at_k")
axC.set_yscale("log"); axC.set_ylim(5e-4, 1.0)
axC.set_xticks(range(len(SYSTEMS))); axC.set_xticklabels(SYSTEMS)
axC.set_xlim(-0.6, len(SYSTEMS) - 0.4)
axC.set_ylabel("Fraction with occupancy ≥ 0.25")
axC.legend(handles=[Line2D([], [], color=COL[n], lw=1.6, label=f"top k by {SLAB[n] if n == 'TCR-IMAP' else 'pair-frequency prior'}") for n, _ in SCORE_KEYS] +
           [Line2D([], [], label="all crystal-negative pairs (base rate)", **REFERENCE_LINE_STYLE)],
           loc="upper center", bbox_to_anchor=(0.5, -0.14), frameon=False, fontsize=5.6, handlelength=1.4,
           ncol=2, columnspacing=1.0)
axC.text(1.0, 0.02, "k = crystal-positive count;  † tie at k", transform=axC.transAxes, ha="right",
         va="bottom", fontsize=5.6, color=TXT)

# ---------------- D: check 1 with docked fractions
axD = fig.add_subplot(gs[1, 1])
xt, xl = [], []
for i, p in enumerate(SYSTEMS):
    for j, r in enumerate(["1", "2", "3", "system"]):
        x = i * 5 + j + (0.3 if r == "system" else 0)
        c1 = chk1[(chk1.pdb == p) & (chk1.rep == r)].iloc[0]
        lw = 2.4 if r == "system" else 1.6
        axD.plot([x - 0.32, x + 0.32], [c1["mean"], c1["mean"]], color="k", lw=lw, solid_capstyle="butt")
        axD.plot(x, c1["median"], marker="x", ms=3.2, mew=0.8, color="k", ls="")
        val("D", f"{p} rep {r} mean / median", f"{c1['mean']:.3f} / {c1['median']:.3f}",
            "md_occupancy_check1_positive.csv", f"pdb={p}, rep={r}")
        xt.append(x); xl.append(r if r != "system" else "all")
        if r != "system":
            dr = dock[(dock.pdb == p) & (dock.rep == int(r))].iloc[0]
            fc, fi = 100 * dr.docked_fraction_count, 100 * dr.docked_fraction_identity
            axD.text(x, 1.02, f"{fc:.0f}", ha="center", va="bottom", fontsize=5.4, color=TXT,
                     transform=axD.get_xaxis_transform())
            axD.text(x, 1.10, f"{fi:.0f}", ha="center", va="bottom", fontsize=5.4, color=TXT,
                     transform=axD.get_xaxis_transform())
            val("D", f"{p} rep {r} docked %", f"{fc:.0f} (count) / {fi:.0f} (identity)",
                "md_occupancy_docked_frames.csv", f"pdb={p}, rep={r}: docked_fraction_count / _identity")
    axD.text(i * 5 + 1.6, -0.12, p, ha="center", va="top", fontsize=7, transform=axD.get_xaxis_transform())
axD.text(-1.2, 1.02, "docked % (count)", ha="right", va="bottom", fontsize=5.4, color=TXT, transform=axD.get_xaxis_transform())
axD.text(-1.2, 1.10, "docked % (identity)", ha="right", va="bottom", fontsize=5.4, color=TXT, transform=axD.get_xaxis_transform())
axD.axhline(0.25, color="0.55", lw=0.6, ls=(0, (3, 2)))
axD.set_xticks(xt); axD.set_xticklabels(xl, fontsize=5.8)
axD.set_ylim(0, 1.0); axD.set_xlim(-1.0, 19.0)
axD.set_ylabel("Crystal-positive occupancy")
axD.legend(handles=[Line2D([], [], color="k", lw=1.6, label="mean"),
                    Line2D([], [], color="k", marker="x", ls="", ms=3.2, mew=0.8, label="median")],
           loc="upper right", frameon=False, fontsize=5.6, handlelength=1.2, borderaxespad=0.1, ncol=2,
           bbox_to_anchor=(1.0, 0.98))

# ---------------- E: backbone RMSD, PBC-corrected characterisation tables
cman = pd.read_csv(T + "md_characterisation_manifest.csv").set_index("item")["value"]
assert cman["occupancy_tables_unchanged"] == "True" and "nojump" in cman["coordinate_measures_trajectory"]
assert cman["md_characterisation_summary.csv_sha256"] == sha(T + "md_characterisation_summary.csv")
assert cman["md_characterisation_timeseries.csv_sha256"] == sha(T + "md_characterisation_timeseries.csv")
csum = pd.read_csv(T + "md_characterisation_summary.csv", **RD)
cts = pd.read_csv(T + "md_characterisation_timeseries.csv", **RD).dropna(subset=["backbone_rmsd"])
axE = []
for ci, p in enumerate(SYSTEMS):
    ax = fig.add_subplot(gsE[0, ci]); axE.append(ax)
    rmeans, fr = [], []
    for r in REP:
        t = cts[(cts.pdb == p) & (cts.rep == r)]
        sm = float(csum[(csum.pdb == p) & (csum.rep == r)].backbone_rmsd_mean.iloc[0])
        assert len(t) == 9001 and np.isclose(t.backbone_rmsd.mean(), sm, rtol=0, atol=1e-9), (p, r)
        ax.plot(t.time_ps / 1000, t.backbone_rmsd, color=REPLICATE_COLORS[r], lw=0.4, rasterized=True, zorder=4 - r)
        rmeans.append(sm)
        fc = 100 * float(dock[(dock.pdb == p) & (dock.rep == r)].docked_fraction_count.iloc[0])
        fr.append(fc)
        val("E", f"{p} rep {r} backbone RMSD mean", f"{sm:.2f} A", "md_characterisation_summary.csv",
            f"pdb={p}, rep={r}: backbone_rmsd_mean (= mean of its 9,001 time-series rows, asserted)")
        val("E", f"{p} rep {r} docked % (count)", f"{fc:.0f}", "md_occupancy_docked_frames.csv",
            f"pdb={p}, rep={r}: docked_fraction_count")
    smean = float(np.mean(rmeans))
    ax.axhline(smean, color="k", lw=0.7, ls=(0, (3, 2)), zorder=5)
    ax.text(197, smean + 0.18, f"mean {smean:.2f} \u00c5", ha="right", va="bottom", fontsize=5.8, color="k", zorder=6,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.6))
    val("E", f"{p} system mean", f"{smean:.2f} A", "md_characterisation_summary.csv",
        f"pdb={p}: mean of the three rep backbone_rmsd_mean cells (equal frame counts)")
    ax.text(0.03, 0.97, "docked % (count), rep 1/2/3: " + " / ".join(f"{v:.0f}" for v in fr), transform=ax.transAxes,
            ha="left", va="top", fontsize=5.5, color=TXT)
    ax.set_xlim(20, 200); ax.set_xticks([50, 100, 150, 200]); ax.set_ylim(0, 8.5)   # the analysed window, 20-200 ns
    ax.set_title(p, fontsize=7, pad=2)
    ax.set_xlabel("Time (ns)", fontsize=7)
    if ci == 0:
        ax.set_ylabel("Backbone RMSD (\u00c5)")
    else:
        ax.tick_params(labelleft=False)
axE[-1].legend(handles=[Line2D([], [], color=REPLICATE_COLORS[r], lw=1.4, label=f"rep {r}") for r in REP] +
               [Line2D([], [], color="k", lw=0.7, ls=(0, (3, 2)), label="system mean")],
               loc="lower right", frameon=False, fontsize=5.6, ncol=2, handlelength=1.3, columnspacing=0.8,
               borderaxespad=0.1)

panel_letters(fig, {"A": [axA], "B": axB_all, "C": [axC], "D": [axD], "E": axE})   # one fixed offset (style.py)

png, pdf = save_fig(fig, OUT_STEM)

# ------------------------------------------------------------------ validation
print(f"\npanel-to-source validation ({len(VAL)} drawn values):")
print_validation_table(VAL, [("panel", "panel", None), ("item", "item", None), ("shown", "drawn", None),
                             ("source", "table", None), ("cell", "cell", None)])
print(f"\nrows flagged not interpretable anywhere in A: {'yes' if any_ni else 'none'}")
n3, arial, bases = pdf_font_check(pdf)
print(f"PDF font check {pdf}: Type3 = {n3}, ArialMT present = {arial}, BaseFonts = {bases}")
assert n3 == 0 and arial
