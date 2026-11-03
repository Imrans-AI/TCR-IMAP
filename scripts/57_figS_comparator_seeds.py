#!/usr/bin/env python3
"""SI figure S6 (p84v2) - the 88-structure comparator subset, per training seed.

The four cells of the comparator-subset analysis (Figure 4F): the comparator's 4.5 A and 6 A cutoffs, each under a
single random split of the kind the comparator used (n = 22 test complexes) and under tcr70 (n = 88). Every model is
an ensemble of three training seeds (predictions averaged, then scored); this figure shows the ensemble and each seed.

A  TCR-IMAP - Struct-TCS per cell: the reported ensemble (cluster-disjoint validation) with 95 % CI and p, the three
   per-seed differences, and the superseded ensemble under the random validation split with its 95 % CI
B  macro-AUPRC of TCR-IMAP and Struct-TCS per cell: ensemble and the three seeds

Sources (nothing is trained, predicted or re-scored):
- results/review/pre_md_run/B_perseed_stability.tsv and B_perseed_macro.tsv (the per-seed table of PRE_MD_RUN_REPORT.md)
- results/tables/_partial/p84v2/31b_structcs_{45A,6A}/{tag}_{setting}.npz: '__percomplex__' (ensemble per-complex AUPRC,
  the arrays behind Figure 4F) and '__perseed__' (per-seed per-complex AUPRC). Every per-seed macro in B_perseed_macro.tsv
  is asserted equal to the mean of its stored list; the ensemble difference is recomputed at 20,000 resamples, seed 42,
  and asserted equal to F_structcs_differences_20k.tsv and to B_perseed_stability.tsv
- results/review/validation_fix/F_structcs_differences_20k.tsv, rows 'random split (superseded)': the superseded
  ensemble at 20,000 resamples (its per-complex arrays were not persisted, so it is drawn from this table and its
  recorded 2,000-resample cell is checked against structcs_randomsplit_deterministic/structcs_setting_p84v2_*.csv)
- results/review/validation_fix/F_structcs_macro.tsv, results/tables/structcs_setting_p84v2_{4p5A,6A}.csv
Fill = 95 % CI excludes zero. Marker = validation split: circle cluster-disjoint (reported), triangle random
(superseded), as in Figure S4. Crosses = single seeds (no test).
"""
import sys, json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, "src")
from viz.style import (apply_style, save_fig, print_validation_table, pstr, sig_marker, signed, minus, panel_letters,
                       PAGE_WIDTH_IN, MODEL_LABEL, METHOD_COLORS)
from viz.figstats import pooled_cluster_bootstrap, ci_excludes_zero, parse_cell, pdf_font_check, N_BOOT, SEED

apply_style()
T, V, P = "results/tables/", "results/review/validation_fix/", "results/review/pre_md_run/"
OUT_STEM = "figures/figS_comparator_subset_seeds_p84v2_v1"
TXT = "0.30"
RD = dict(float_precision="round_trip")
SEEDS = (42, 43, 44)
CELLS = [("4p5A", "their_random"), ("4p5A", "our_tcr70"), ("6A", "their_random"), ("6A", "our_tcr70")]
CUT = {"4p5A": "4.5 Å", "6A": "6 Å"}
SET = {"their_random": "single random split", "our_tcr70": "tcr70"}
CACHE = {"4p5A": "31b_structcs_45A", "6A": "31b_structcs_6A"}
VAL = []
def val(panel, what, shown, src, cell):
    VAL.append(dict(panel=panel, what=what, shown=shown, src=src, cell=cell))

stab = pd.read_csv(P + "B_perseed_stability.tsv", sep="\t", **RD).set_index(["cutoff", "setting"])
pmac = pd.read_csv(P + "B_perseed_macro.tsv", sep="\t", **RD)
F = pd.read_csv(V + "F_structcs_differences_20k.tsv", sep="\t", **RD)
FM = pd.read_csv(V + "F_structcs_macro.tsv", sep="\t", **RD)
D = {}
for tag, st in CELLS:
    z = np.load(T + f"_partial/p84v2/{CACHE[tag]}/{tag}_{st}.npz", allow_pickle=True)
    pc, ps = json.loads(str(z["__percomplex__"])), json.loads(str(z["__perseed__"]))
    o, s, cl = np.array(pc["Ours_auprc"]), np.array(pc["Struct-TCS_auprc"]), np.array(pc["clu70"])
    ens = pooled_cluster_bootstrap(o, s, cl, N_BOOT, SEED)
    srow = stab.loc[(tag, st)]
    frow = F[(F.cutoff == tag) & (F.setting == st) & (F.comparison == "TCR-IMAP vs Struct-TCS")].set_index("config")
    rep, sup = frow.loc["cluster-disjoint (reported)"], frow.loc["random split (superseded)"]
    for f in ("diff", "lo", "hi", "p"):
        assert abs(ens[f] - rep[f]) < 1e-12, (tag, st, f)
    assert abs(ens["diff"] - srow.ensemble) < 1e-12 and abs(ens["p"] - srow.ens_p) < 1e-12 and int(srow.n_test) == len(o)
    # superseded ensemble: 20k value from the table; its recorded 2k cell must equal the committed random-split csv
    rs = pd.read_csv(V + f"structcs_randomsplit_deterministic/structcs_setting_p84v2_{tag}.csv").set_index("setting")
    assert parse_cell(sup.recorded_2k) == parse_cell(rs.loc[st, "vs_Struct-TCS"]), (tag, st)
    # per seed: macro = mean of the stored per-complex list; difference = the stability table's
    seeds = {}
    for sd in SEEDS:
        mo, ms = (float(np.mean(ps[f"{m}_seed{sd}_auprc"])) for m in ("Ours", "Struct-TCS"))
        r = pmac[(pmac.cutoff == tag) & (pmac.setting == st) & (pmac.seed == sd)].iloc[0]
        assert abs(mo - r.TCR_IMAP) < 1e-12 and abs(ms - r.Struct_TCS) < 1e-12, (tag, st, sd)
        assert abs((mo - ms) - srow[f"seed{sd}"]) < 1e-12
        seeds[sd] = (mo, ms, mo - ms)
    # ensemble macro means
    mrow = FM[(FM.cutoff == tag) & (FM.setting == st) & (FM.config == "cluster-disjoint (reported)")].set_index("method")
    assert abs(np.nanmean(o) - mrow.loc["TCR-IMAP", "macro_auprc"]) < 1e-12 and abs(np.nanmean(s) - mrow.loc["Struct-TCS", "macro_auprc"]) < 1e-12
    signs = {np.sign(v[2]) for v in seeds.values()}
    assert len(signs) == 1 and bool(srow.sign_flip) is False, (tag, st, "seeds disagree in sign")
    D[(tag, st)] = dict(ens=ens, sup=dict(diff=sup["diff"], lo=sup.lo, hi=sup.hi, p=sup.p), seeds=seeds,
                        macro=(float(np.nanmean(o)), float(np.nanmean(s))), n=len(o))
r6 = D[("6A", "their_random")]
assert r6["sup"]["diff"] > 0 and r6["ens"]["diff"] < 0 and r6["ens"]["hi"] < 0, "the 6 A single-random-split reversal"
print(f"4 cells: ensemble differences recomputed at {N_BOOT:,} resamples, seed {SEED} = F_structcs_differences_20k.tsv and "
      "B_perseed_stability.tsv; 12 per-seed macros = means of the stored per-seed lists; superseded 2k cells = committed "
      "random-split csv; three seeds agree in sign in all 4 cells; 6 A single random split reverses (+ -> -, CI below 0)")

# ------------------------------------------------------------------ figure
FW, FH = PAGE_WIDTH_IN, 3.70
fig = plt.figure(figsize=(FW, FH))
cS, cM = METHOD_COLORS["Struct-TCS"], METHOD_COLORS[MODEL_LABEL]
yy = np.array([3.0, 2.0, 0.6, -0.4])        # 4.5 A cells, gap, 6 A cells
ylab = [f"{CUT[t]}, {SET[s]}\nn = {D[(t, s)]['n']}" for t, s in CELLS]

# --- A
LA, WA = 1.45, 1.70
ax = axA = fig.add_axes([LA / FW, 1.00 / FH, WA / FW, 2.25 / FH])
XD, XP, XS = 1 + 0.10 / WA, 1 + 1.40 / WA, 1 + 1.98 / WA
for (tag, st), y in zip(CELLS, yy):
    d = D[(tag, st)]
    for dy, r, mk in ((+0.16, d["ens"], "o"), (-0.18, d["sup"], "^")):
        ax.plot([r["lo"], r["hi"]], [y + dy, y + dy], color=cS, lw=1.0 if mk == "o" else 0.8,
                solid_capstyle="butt", alpha=1.0 if mk == "o" else 0.55)
        ax.plot(r["diff"], y + dy, mk, ms=4.6 if mk == "o" else 4.0, zorder=5, alpha=1.0 if mk == "o" else 0.55,
                **sig_marker(cS, ci_excludes_zero(r["lo"], r["hi"])))
    for sd, (_, _, dd) in d["seeds"].items():
        ax.plot(dd, y + 0.16, "x", ms=3.8, mew=0.8, color="0.25", zorder=4)
    e = d["ens"]
    ax.text(XD, y + 0.16, minus(f"{signed(e['diff'], 4)} [{signed(e['lo'], 4)}, {signed(e['hi'], 4)}]"),
            transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.7, color=TXT)
    ax.text(XP, y + 0.16, pstr(e["p"]), transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.7, color=TXT)
    ax.text(XS, y + 0.16, minus(", ".join(signed(v[2], 3) for v in d["seeds"].values())),
            transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.7, color=TXT)
    sp = d["sup"]
    ax.text(XD, y - 0.18, minus(f"{signed(sp['diff'], 4)} [{signed(sp['lo'], 4)}, {signed(sp['hi'], 4)}]"),
            transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.5, color="0.55")
    ax.text(XP, y - 0.18, pstr(sp["p"]), transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.5, color="0.55")
    val("A", f"{CUT[tag]} {SET[st]} ensemble (reported)", f"{signed(e['diff'], 4)} [{signed(e['lo'], 4)}, {signed(e['hi'], 4)}] "
        f"{pstr(e['p'])} ({'filled' if ci_excludes_zero(e['lo'], e['hi']) else 'hollow'})", "F_structcs_differences_20k.tsv",
        f"cutoff={tag}, setting={st}, Struct-TCS, cluster-disjoint (reported)")
    val("A", f"{CUT[tag]} {SET[st]} seeds 42/43/44", ", ".join(f"{v[2]:+.4f}" for v in d["seeds"].values()),
        "B_perseed_stability.tsv", f"cutoff={tag}, setting={st}, seed42/43/44")
    val("A", f"{CUT[tag]} {SET[st]} ensemble, random validation (superseded)",
        f"{signed(sp['diff'], 4)} [{signed(sp['lo'], 4)}, {signed(sp['hi'], 4)}] {pstr(sp['p'])}",
        "F_structcs_differences_20k.tsv", f"cutoff={tag}, setting={st}, Struct-TCS, random split (superseded)")
ax.axvline(0, ls="--", lw=0.8, color="0.4")
ax.set_yticks(yy); ax.set_yticklabels(ylab, fontsize=6.2, linespacing=1.15)
ax.set_ylim(-0.95, 3.55)
ax.set_xlabel(f"{MODEL_LABEL} − Struct-TCS\n(per-residue macro-AUPRC)", fontsize=7, linespacing=1.2)
for x, h in ((XD, "Δ [95% CI]"), (XP, "p"), (XS, "seeds 42, 43, 44")):
    ax.text(x, 1.0, h, transform=ax.transAxes, ha="left", va="bottom", fontsize=6.2, color="0.2")
fig.legend(handles=[Line2D([], [], marker="o", ls="", mfc=cS, mec=cS, ms=4.3, label="ensemble, cluster-disjoint validation (reported)"),
                   Line2D([], [], marker="x", ls="", color="0.25", ms=3.8, mew=0.8, label="single seed"),
                   Line2D([], [], marker="^", ls="", mfc=cS, mec=cS, ms=4.0, alpha=0.55,
                          label="ensemble, random validation (superseded)")],
          loc="lower left", bbox_to_anchor=(0.01, 0.06), ncol=3, frameon=False, fontsize=6.2, handletextpad=0.3,
          columnspacing=1.0)
ax.set_title("Paired difference per cell", fontsize=8, pad=12)

# --- B
ax = axB = fig.add_axes([6.08 / FW, 1.00 / FH, 0.95 / FW, 2.25 / FH])
for (tag, st), y in zip(CELLS, yy):
    d = D[(tag, st)]
    for (mi, col, lab), dy in (((0, cM, MODEL_LABEL), +0.14), ((1, cS, "Struct-TCS"), -0.14)):
        ax.plot(d["macro"][mi], y + dy, "o", ms=4.0, color=col, mec="white", mew=0.4, zorder=3)
        for sd, v in d["seeds"].items():
            ax.plot(v[mi], y + dy, "x", ms=3.6, mew=0.8, color=col, zorder=2)
        val("B", f"{CUT[tag]} {SET[st]} {lab}", f"ensemble {d['macro'][mi]:.4f}; seeds " +
            ", ".join(f"{v[mi]:.4f}" for v in d["seeds"].values()), "F_structcs_macro.tsv; B_perseed_macro.tsv",
            f"cutoff={tag}, setting={st}, {lab}")
ax.set_yticks(yy); ax.set_yticklabels([])
ax.set_ylim(-0.95, 3.55); ax.set_xlim(0.55, 0.95)
ax.set_xlabel("Macro-AUPRC", fontsize=7)
fig.legend(handles=[Line2D([], [], marker="o", ls="", color=c, mec="white", ms=4.2, label=l)
                    for l, c in ((MODEL_LABEL, cM), ("Struct-TCS", cS))] +
           [Line2D([], [], marker="x", ls="", color="0.35", ms=3.6, mew=0.8, label="single seed")],
           loc="lower left", bbox_to_anchor=(5.80 / FW, 0.06), ncol=1, frameon=False, fontsize=6.2, handletextpad=0.2)
ax.set_title("Per method", fontsize=8, pad=12)

fig.text(0.99, 0.008, minus(f"A: pooled clu70 cluster bootstrap, {N_BOOT:,} resamples, seed {SEED}; the reported ensemble "
                            "recomputed from stored per-complex arrays"), ha="right", va="bottom", fontsize=5.4, color="0.45")
panel_letters(fig, {"A": [axA], "B": [axB]})
png, pdf = save_fig(fig, OUT_STEM)
plt.close(fig)

print("\n--- every drawn number and its source ---")
print_validation_table(VAL, [("panel", "panel", None), ("what", "value", None), ("shown", "drawn", None),
                             ("src", "source", None), ("cell", "cell", None)])
n3, arial, bases = pdf_font_check(pdf)
print(f"\nPDF font check {pdf}: Type3 = {n3}, ArialMT present = {arial}, BaseFonts = {bases}")
assert n3 == 0 and arial
