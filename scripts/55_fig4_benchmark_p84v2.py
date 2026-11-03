#!/usr/bin/env python3
"""Figure 4 (p84v2) - benchmark against published comparators, tcr70 primary.

Built from 37d_fig4_benchmark_p82v1.py; inputs are the p84v2 tables (insertion-code-keyed labels,
position-aligned tokens; results/review/CORRECTION_PLAN_P84V2.md). Under p84v2 TCR-IMAP is ahead of
Struct-TCS in every protocol; random (+0.0066) and temporal (+0.021) are annotated as not significant.

A  per-residue macro-AUPRC by protocol (random, double, allele n = 64, temporal, tcr70) for
   TCR-IMAP, Struct-TCS, Seq-TCS, position prior and NN transfer. The protocols where the
   difference to Struct-TCS is not significant (random, temporal) are annotated, not hidden.
B  paired difference TCR-IMAP - Struct-TCS per protocol, 95 % CI
C  pairwise tcr70, per complex: TCR-IMAP vs pair-frequency prior
D  pairwise tcr70, per complex: TCR-IMAP vs TEIM-Res
E  cluster-macro sensitivity (69 clu70 clusters equally weighted) next to complex-macro
F  88-structure comparator subset at 4.5 A and 6 A, 20,000-resample clu70 cluster bootstrap from the 31b per-complex
   arrays -- shared subset and cutoff only; this is NOT the comparator's own evaluation (candidate positions, TCR representation, pooled task
   and validation all differ)

Statistics. Per-complex macro-AUPRC. Wherever per-complex arrays exist (A-E) the CI and p are
recomputed here by the pooled clu70 cluster bootstrap at 20,000 resamples, seed 42, and those
are displayed; the recorded 2,000-resample values are first reproduced exactly (seed 0) as a
provenance check, and both are printed in the validation table. TEIM-Res has no recorded
difference (36b prints only its macro-AUPRC), so D is recomputed without a 2k comparison. Cluster-macro (E): cluster_macro_p84v2.csv
was written by scripts/55_stats_p84v2.py with the estimator verified on p82v1; it is reproduced
here at 2,000 resamples, seed 0, before the 20k recompute.
F uses the per-complex arrays the 31b correction persisted (results/tables/_partial/p84v2/31b_structcs_*), at
20,000 resamples, seed 42; their macro means reproduce the recorded table, which is printed alongside.
p is displayed only through style.pstr(). No Wilcoxon/Friedman across complexes.

Encoding. Difference markers take the COMPARATOR's colour (METHOD_COLORS); fill = 95 % CI
excludes zero (sig_marker), hollow = CI includes zero and nothing else.
"""
import sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, "."); sys.path.insert(0, "src")
from viz.style import (apply_style, save_fig, print_validation_table, pstr, sig_marker, signed, minus,
                       PAGE_WIDTH_IN, PANEL_LABEL_STYLE as PL, MODEL_LABEL, METHOD_COLORS, panel_letters)
from viz.figstats import (pooled_cluster_bootstrap as pcb, cluster_macro_bootstrap as cmb,
                          parse_cell, ci_excludes_zero, pdf_font_check, N_BOOT, SEED)

apply_style()
T = "results/tables/"
OUT_STEM = "figures/fig4_benchmark_p84v2_v1"
PROTOS = ["random", "double", "allele", "temporal", "tcr70"]       # keys as stored in the frozen files
# "double" is a legacy intersection-selected partition, not a double holdout: its held-out groups
# remain in training (results/review/V3_PART_A_VERIFICATION.md A1). The protocols are NOT ordered
# by strictness; they probe different distribution shifts on different subsets.
PLAB = {"random": "random", "double": "intersection\n(legacy)", "allele": "allele",
        "temporal": "temporal", "tcr70": "tcr70"}
KEY2LAB = {"Ours": MODEL_LABEL, "Struct-TCS": "Struct-TCS", "Seq-TCS": "Seq-TCS",
           "PosFreq": "Position prior", "NN": "NN transfer"}
METHODS = list(KEY2LAB)
C = {k: METHOD_COLORS[v] for k, v in KEY2LAB.items()}
C["pair"] = METHOD_COLORS["Pair-frequency prior"]; C["TEIM"] = METHOD_COLORS["TEIM-Res"]
TXT = "0.30"
FLOOR2K = 1 / 2000

# ------------------------------------------------------------------ data + provenance
bench = pd.read_csv(T + "benchmark_p84v2.csv").set_index("protocol")
pcx = {p: np.load(T + f"percomplex_{p}_p84v2.npz", allow_pickle=True) for p in PROTOS}
for p in PROTOS:
    assert len(pcx[p]["pdb"]) == int(bench.loc[p, "n"]), p
    for k in METHODS:
        assert round(float(np.nanmean(pcx[p][f"{k}_auprc"])), 4) == bench.loc[p, k], (p, k)
assert len(pcx["allele"]["pdb"]) == 64

po = np.load(T + "pairwise_ours_p84v2.npz", allow_pickle=True)
te = np.load(T + "teim_res_pairwise_p84v2.npz", allow_pickle=True)
pr = pd.read_csv(T + "percomplex_tcr70_paired_p84v2.tsv", sep="\t")
assert (po["pdb"] == te["pdb"]).all() and (po["pdb"] == pr.pdb.values).all()
assert np.allclose(po["auprc"], pr.grip.values) and (po["clu70"] == pr.clu70.values).all()
md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
assert (md.loc[po["pdb"], "clu70"].values == po["clu70"]).all() and md.clu70.nunique() == 69
pw_model, pw_prior, pw_teim, pw_cl = po["auprc"], pr.prior.values, te["auprc"], po["clu70"]

cm = pd.read_csv(T + "cluster_macro_p84v2.csv")
sub = {tag: pd.read_csv(T + f"structcs_setting_p84v2_{tag}.csv").set_index("setting")
       for tag in ("4p5A", "6A")}
# Panel F at the paper's standard (20,000 resamples, seed 42). The 31b runner (53_structcs_setting_p84v2.py) persisted
# the per-complex AUPRC arrays it used for the recorded table; their macro means must reproduce that table, and the
# 20k result must equal the values reported in FINISH_PENDING_REPORT.md (31b addendum).
import json as _json
F20 = {}
for tag, cache in (("4p5A", "31b_structcs_45A"), ("6A", "31b_structcs_6A")):
    for st_ in ("their_random", "our_tcr70"):
        _z = np.load(T + f"_partial/p84v2/{cache}/{tag}_{st_}.npz", allow_pickle=True)
        _pc = _json.loads(str(_z["__percomplex__"]))
        _o, _s, _c = np.array(_pc["Ours_auprc"]), np.array(_pc["Struct-TCS_auprc"]), np.array(_pc["clu70"])
        _row = sub[tag].loc[st_]
        assert len(_o) == int(_row.n_test), (tag, st_)
        assert round(float(np.nanmean(_o)), 4) == float(_row.Ours_auprc) and \
               round(float(np.nanmean(_s)), 4) == float(_row["Struct-TCS_auprc"]), (tag, st_, "macro means differ")
        F20[(tag, st_)] = pcb(_o, _s, _c)
F_REPORTED = {("4p5A", "their_random"): (+0.0269, -0.0102, +0.0707), ("4p5A", "our_tcr70"): (+0.1057, +0.0705, +0.1357),
              ("6A", "their_random"): (-0.0949, -0.1335, -0.0517), ("6A", "our_tcr70"): (+0.0964, +0.0684, +0.1235)}
for k, (d_, lo_, hi_) in F_REPORTED.items():
    r_ = F20[k]
    assert (round(r_["diff"], 4), round(r_["lo"], 4), round(r_["hi"], 4)) == (d_, lo_, hi_), (k, r_)

VAL = []   # validation rows


def ptab(p, n):
    """p for the validation table: never print a floored value as a number."""
    return f"floor ({pstr(p)})" if p <= 1 / n + 1e-12 else f"{p:.4f}"


def compare(a, b, cl, rec, panel, name, macro=False):
    """Reproduce the recorded 2,000-resample result (seed 0), then recompute at 20k/seed 42."""
    f = cmb if macro else pcb
    r2, r20 = f(a, b, cl, 2000, 0), f(a, b, cl, N_BOOT, SEED)
    if rec is None:      # nothing recorded to reproduce: show the seed-0 2k value for reference
        rep, rec = None, (r2["diff"], r2["lo"], r2["hi"], r2["p"])
    else:
        rep = (all(abs(round(x, 4) - y) < 1e-9 for x, y in zip((r2["diff"], r2["lo"], r2["hi"]), rec[:3]))
               and abs(round(r2["p"], 4) - rec[3]) < 1e-9)
    VAL.append(dict(panel=panel, comparison=name, n=f"{r20['n']}/{r20['k']}",
                    recorded=f"{rec[0]:+.4f} [{rec[1]:+.4f},{rec[2]:+.4f}] {ptab(rec[3], 2000)}",
                    re2k="none recorded" if rep is None else ("exact" if rep else "MISMATCH"),
                    r20=f"{r20['diff']:+.4f} [{r20['lo']:+.4f},{r20['hi']:+.4f}] {ptab(r20['p'], N_BOOT)}",
                    shown=pstr(r20["p"]), fill="filled" if ci_excludes_zero(r20["lo"], r20["hi"]) else "hollow"))
    return r20


# ------------------------------------------------------------------ statistics
perres = {}
for p in PROTOS:
    z = pcx[p]
    for k in ("Struct-TCS", "Seq-TCS", "PosFreq", "NN"):
        perres[(p, k)] = compare(z["Ours_auprc"], z[f"{k}_auprc"], z["clu70"],
                                 parse_cell(bench.loc[p, f"vs_{k}"]),
                                 "A/B" if k == "Struct-TCS" else "A", f"per-res {p}: vs {KEY2LAB[k]}")
cmr = cm.set_index("model")
# recorded complex-macro CI/p for the pairwise comparison live in protocol_main_p84v2.csv
_pm = pd.read_csv(T + "protocol_main_p84v2.csv").set_index("protocol").loc["tcr70"]
pw_vs_prior = compare(pw_model, pw_prior, pw_cl,
                      (round(float(_pm["diff"]), 4), float(_pm.ci_lo), float(_pm.ci_hi), float(_pm.pval)),
                      "C/E", "pairwise tcr70: vs pair-frequency prior")
pw_vs_teim = compare(pw_model, pw_teim, pw_cl, None, "D", "pairwise tcr70: vs TEIM-Res")
# p84v2: the 36b runner logs only TEIM-Res macro-AUPRC (results/review/train_p84v2_logs/teim_res.log);
# no recorded difference exists, so there is nothing to reproduce

z70 = pcx["tcr70"]
macro = {}
for key, (a, b, cl, mdl) in {
        "Struct-TCS": (z70["Ours_auprc"], z70["Struct-TCS_auprc"], z70["clu70"], "Ours_minus_Struct-TCS"),
        "PosFreq": (z70["Ours_auprc"], z70["PosFreq_auprc"], z70["clu70"], "Ours_minus_PosFreq"),
        "pair": (pw_model, pw_prior, pw_cl, "model_minus_freq")}.items():
    r = cmr.loc[mdl]
    rec = (r.paired_diff_vs_GRIP_cluster_macro, r.ci_lo_cluster_macro, r.ci_hi_cluster_macro, r.p_cluster_macro)
    macro[key] = compare(a, b, cl, rec, "E", f"cluster-macro tcr70: vs {key}", macro=True)
macro_cx = {"Struct-TCS": perres[("tcr70", "Struct-TCS")], "PosFreq": perres[("tcr70", "PosFreq")],
            "pair": pw_vs_prior}

win_prior = int((pw_model > pw_prior).sum()); win_teim = int((pw_model > pw_teim).sum())
n_pw = len(pw_model)

# ------------------------------------------------------------------ figure
# explicit axes positions in inches: constrained layout collapses with equal-aspect scatters
# and p-value text placed outside the forest axes
FW, FH = PAGE_WIDTH_IN, 7.85
fig = plt.figure(figsize=(FW, FH))
def axin(l, b, w, h): return fig.add_axes([l / FW, b / FH, w / FW, h / FH])
PANELS = {}
def letter(ax, lab, dx_in=None):     # registers the panel; style.panel_letters places every letter at one offset
    PANELS.setdefault(lab, []).append(ax)

# --- A ---------------------------------------------------------------------------------
ax = axin(0.55, 5.65, 4.05, 1.85)
OFF = np.linspace(-0.30, 0.30, len(METHODS))
for i, p in enumerate(PROTOS):
    for o, k in zip(OFF, METHODS):
        v = float(np.nanmean(pcx[p][f"{k}_auprc"]))
        ax.plot(i + o, v, "o", ms=5 if k == "Ours" else 4.2, color=C[k], mec="white", mew=0.4,
                zorder=4 if k == "Ours" else 3)
    if i < len(PROTOS) - 1:
        ax.axvline(i + 0.5, color="0.88", lw=0.6, zorder=0)
ax.set_xticks(range(len(PROTOS)))
ax.set_xticklabels([f"{PLAB[p]}\nn = {len(pcx[p]['pdb'])}" + ("\n(primary)" if p == "tcr70" else "")
                    for p in PROTOS], fontsize=6.5, linespacing=1.25)
ax.set_xlim(-0.5, len(PROTOS) - 0.5); ax.set_ylim(0.30, 1.03)
ax.set_ylabel("Per-residue macro-AUPRC")
NS_PROTOS = [p for p in PROTOS if not ci_excludes_zero(perres[(p, "Struct-TCS")]["lo"], perres[(p, "Struct-TCS")]["hi"])]
# Under the cluster-disjoint validation split (CORRECTION_PLAN_P84V2.md, "Per-residue validation split"),
# tcr70 joins random and temporal: its 95 % CI now includes zero (+0.0312, p = 0.088).
assert NS_PROTOS == ["random", "temporal", "tcr70"], NS_PROTOS
# Every protocol cell is annotated the same way, with the two plotted means only. The paired test, and so the
# significance verdict, is panel B's for all five protocols (difference, p, fill), so A carries no verdict text.
for i, p in enumerate(PROTOS):
    mo, ms_ = float(np.nanmean(pcx[p]["Ours_auprc"])), float(np.nanmean(pcx[p]["Struct-TCS_auprc"]))
    ax.text(i, 1.02, minus(f"{MODEL_LABEL} {mo:.3f}\nStruct-TCS {ms_:.3f}"), ha="center", va="top",
            fontsize=5.6, color=TXT, linespacing=1.25)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=C[k], mec="white", mew=0.4, ms=4.5,
                          label=KEY2LAB[k]) for k in METHODS],
          loc="lower left", ncol=5, frameon=False, fontsize=6, handletextpad=0.2, columnspacing=0.9,
          bbox_to_anchor=(0.0, -0.01))
ax.set_title("Per-residue benchmark across protocols", fontsize=8, pad=3)
letter(ax, "A", 0.50)

# --- B ---------------------------------------------------------------------------------
ax = axin(5.25, 5.65, 1.85, 1.85)
cS = C["Struct-TCS"]
for i, p in enumerate(PROTOS):
    r = perres[(p, "Struct-TCS")]
    ax.plot([i, i], [r["lo"], r["hi"]], color=cS, lw=1.0, solid_capstyle="butt")
    ax.plot(i, r["diff"], "o", ms=4.8, zorder=3, **sig_marker(cS, ci_excludes_zero(r["lo"], r["hi"])))
    ax.text(i, r["hi"] + 0.004, f"{signed(r['diff'], 4)}\n{pstr(r['p'])}", ha="center", va="bottom",
            fontsize=5.4, color=TXT, linespacing=1.2)
ax.axhline(0, ls="--", lw=0.8, color="0.4")
ax.set_xticks(range(len(PROTOS))); ax.set_xticklabels([PLAB[p].replace("\n", " ") for p in PROTOS], fontsize=6.0, rotation=35, ha="right",
                                                       rotation_mode="anchor")
_bl = [perres[(p, "Struct-TCS")] for p in PROTOS]
ax.set_xlim(-0.5, len(PROTOS) - 0.5); ax.set_ylim(min(-0.065, min(r["lo"] for r in _bl) - 0.01),
                                               max(0.105, max(r["hi"] for r in _bl) + 0.035))
ax.set_ylabel(f"{MODEL_LABEL} − Struct-TCS\n(per-residue macro-AUPRC)", linespacing=1.2)
ax.set_title("Paired difference vs Struct-TCS", fontsize=8, pad=3)
letter(ax, "B", 0.62)

# --- C / D -----------------------------------------------------------------------------
def scatter(ax, comp, ccol, clab, res, wins, lab_):
    lim = 0.92
    ax.plot([0, lim], [0, lim], ls="--", lw=0.7, color="0.6", zorder=1)
    ax.scatter(comp, pw_model, s=9, color=ccol, alpha=0.75, edgecolor="none", zorder=2)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
    ax.set_xlabel(f"{clab[0].upper() + clab[1:]} pairwise AUPRC"); ax.set_ylabel(f"{MODEL_LABEL} pairwise AUPRC")
    ax.text(0.03, 0.97, f"{MODEL_LABEL} higher on {wins}/{n_pw} ({100 * wins / n_pw:.1f}%)\n"
            f"pooled Δ = {signed(res['diff'], 4)}\n95% CI [{signed(res['lo'], 4)}, {signed(res['hi'], 4)}]\n"
            f"{pstr(res['p'])}", transform=ax.transAxes, ha="left", va="top", fontsize=6,
            color=TXT, linespacing=1.35)
    ax.set_title(f"Pairwise, tcr70: vs {clab}", fontsize=8, pad=3)
    letter(ax, lab_, 0.62)

scatter(axin(0.75, 2.50, 2.40, 2.40), pw_prior, C["pair"], "pair-frequency prior", pw_vs_prior, win_prior, "C")
scatter(axin(4.45, 2.50, 2.40, 2.40), pw_teim, C["TEIM"], "TEIM-Res", pw_vs_teim, win_teim, "D")

# --- E ---------------------------------------------------------------------------------
ax = axin(1.15, 0.45, 1.95, 1.40)
ROWS = [("Struct-TCS", "per-residue\nvs Struct-TCS", C["Struct-TCS"]),
        ("PosFreq", "per-residue\nvs position prior", C["PosFreq"]),
        ("pair", "pairwise\nvs pair-frequency prior", C["pair"])]
yy = np.arange(len(ROWS))[::-1]
for (key, lab, col), y in zip(ROWS, yy):
    for dy, res, mk in ((+0.16, macro_cx[key], "o"), (-0.16, macro[key], "s")):
        ax.plot([res["lo"], res["hi"]], [y + dy, y + dy], color=col, lw=1.0, solid_capstyle="butt")
        ax.plot(res["diff"], y + dy, mk, ms=4.3, zorder=3, **sig_marker(col, ci_excludes_zero(res["lo"], res["hi"])))
        ax.text(1.0, y + dy, pstr(res["p"]), transform=ax.get_yaxis_transform(), ha="left",
                va="center", fontsize=5.6, color=TXT)
ax.axvline(0, ls="--", lw=0.8, color="0.4")
ax.set_yticks(yy); ax.set_yticklabels([r[1] for r in ROWS], fontsize=6.2, linespacing=1.15)
ax.set_ylim(-0.6, len(ROWS) + 0.15)
ax.set_xlabel(f"{MODEL_LABEL} − comparator (macro-AUPRC)")
ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc="0.3", mec="0.3", ms=4, label="complex-macro (174)"),
                   Line2D([], [], marker="s", ls="", mfc="0.3", mec="0.3", ms=4, label="cluster-macro (69)")],
          loc="upper right", ncol=2, frameon=False, fontsize=5.6, handletextpad=0.2,
          columnspacing=0.8, borderaxespad=0.1)
ax.set_title("Cluster-macro sensitivity, tcr70", fontsize=8, pad=3)
letter(ax, "E", 1.10)

# --- F ---------------------------------------------------------------------------------
ax = axin(4.90, 0.45, 1.70, 1.40)
# row labels per author decision 9; the cutoff is a group header inside the axes
FROWS = [("4p5A", "their_random", "single random split\n(comparator design)"),
         ("4p5A", "our_tcr70", "tcr70"),
         ("6A", "their_random", "single random split\n(comparator design)"),
         ("6A", "our_tcr70", "tcr70")]
yy = np.array([3.4, 2.4, 1.0, 0.0])
fvals = []
for (tag, st, lab), y in zip(FROWS, yy):
    row = sub[tag].loc[st]
    _r20 = F20[(tag, st)]
    d, lo, hi, p = _r20["diff"], _r20["lo"], _r20["hi"], _r20["p"]
    fvals.append((tag, st, int(row.n_test), d, lo, hi, p, int(row.n_shared_with_comparator)))
    ax.plot([lo, hi], [y, y], color=cS, lw=1.0, solid_capstyle="butt")
    ax.plot(d, y, "o", ms=4.5, zorder=3, **sig_marker(cS, ci_excludes_zero(lo, hi)))
    ax.text(1.0, y, pstr(p), transform=ax.get_yaxis_transform(), ha="left", va="center",
            fontsize=5.6, color=TXT)
ax.axvline(0, ls="--", lw=0.8, color="0.4")
ax.set_yticks(yy)
ax.set_yticklabels([f"{lab}, n = {v[2]}" for (_, _, lab), v in zip(FROWS, fvals)], fontsize=6.2,
                   linespacing=1.15)
ax.set_ylim(-1.0, 4.3)
ax.set_xlim(min(-0.012, min(v[4] for v in fvals) - 0.008), max(0.115, max(v[5] for v in fvals) + 0.008))
for yh, cut in ((3.95, "4.5 \u00c5 cutoff"), (1.55, "6 \u00c5 cutoff")):
    ax.text(-0.04, yh, cut, transform=ax.get_yaxis_transform(), ha="right", va="center",
            fontsize=6.4, color="0.1")
ax.set_xlabel(f"{MODEL_LABEL} − Struct-TCS (per-residue macro-AUPRC)")
ax.set_title("88 comparator structures:\nshared subset and cutoff only", fontsize=7.5, pad=3,
             linespacing=1.15)
ax.text(1.0, -0.72, f"clu70 cluster bootstrap, {N_BOOT:,} resamples", transform=ax.get_yaxis_transform(),
        ha="right", va="center", fontsize=5.4, color="0.45")
letter(ax, "F", 1.20)

panel_letters(fig, PANELS)          # one fixed offset for every figure (style.py)
png, pdf = save_fig(fig, OUT_STEM)
plt.close(fig)

# ------------------------------------------------------------------ validation
print("\n--- comparisons: recorded (2,000 resamples, reproduced at seed 0) vs recomputed "
      f"({N_BOOT:,} resamples, seed {SEED}); n = complexes/clusters ---")
print_validation_table(VAL, [("panel", "panel", None), ("comparison", "comparison", None),
                             ("n", "n/k", None), ("recorded", "recorded 2k", None),
                             ("re2k", "2k repro", None), ("r20", "recomputed 20k", None),
                             ("shown", "displayed p", None), ("fill", "marker", None)])
print("\n--- panel A: per-residue macro-AUPRC (percomplex_*_p84v2.npz == benchmark_p84v2.csv) ---")
print(bench[["n"] + METHODS].rename(columns=KEY2LAB).to_string())
print(f"\n--- panel F: 31b per-complex arrays, {N_BOOT:,} resamples, seed {SEED} (macro means reproduce the recorded "
      "table; values equal the 31b addendum) ---")
for tag, st, n, d, lo, hi, p, nsh in fvals:
    print(f"  {tag:4s} {st:12s} n_test={n:3d} shared={nsh}  vs Struct-TCS {d:+.4f} [{lo:+.4f},{hi:+.4f}] "
          f"{ptab(p, N_BOOT)} -> shown {pstr(p)} | recorded 2k: {sub[tag].loc[st, 'vs_Struct-TCS']}")
print(f"\nwin rates: vs pair prior {win_prior}/{n_pw} ({100*win_prior/n_pw:.1f}%), "
      f"vs TEIM-Res {win_teim}/{n_pw} ({100*win_teim/n_pw:.1f}%)")
print(f"pairwise macro: TCR-IMAP {np.nanmean(pw_model):.4f}, prior {np.nanmean(pw_prior):.4f}, "
      f"TEIM-Res {np.nanmean(pw_teim):.4f}")
for k in macro:
    print(f"cluster-macro {k}: TCR-IMAP {macro[k]['macro_a']:.4f} vs {macro[k]['macro_b']:.4f}")

# cross-check against the step-1 table (scripts/55_stats_p84v2.py, 20,000 resamples, seed 42)
S20 = pd.read_csv("results/review/figures_p84v2/stats_20k_seed42.tsv", sep="\t").set_index(["group", "comparison"])
def same(res, grp, comp):
    r = S20.loc[(grp, comp)]
    return all(abs(res[k] - r[k]) < 1e-12 for k in ("diff", "lo", "hi", "p"))
chk = [(f"per-res {p}: vs {KEY2LAB[k]}", same(perres[(p, k)], "per-residue", f"{p}: TCR-IMAP vs {k}"))
       for p in PROTOS for k in ("Struct-TCS", "Seq-TCS", "PosFreq", "NN")]
chk += [("pairwise tcr70 vs pair prior", same(pw_vs_prior, "pairwise", "tcr70: TCR-IMAP vs pair-frequency prior")),
        ("pairwise tcr70 vs TEIM-Res", same(pw_vs_teim, "pairwise", "tcr70: TCR-IMAP vs TEIM-Res")),
        ("cluster-macro vs Struct-TCS", same(macro["Struct-TCS"], "cluster-macro", "per-residue tcr70: TCR-IMAP vs Struct-TCS")),
        ("cluster-macro vs position prior", same(macro["PosFreq"], "cluster-macro", "per-residue tcr70: TCR-IMAP vs PosFreq")),
        ("cluster-macro vs pair prior", same(macro["pair"], "cluster-macro", "pairwise tcr70: TCR-IMAP vs pair-frequency prior")),
        ("win rate vs prior == step-1", f"{win_prior}/{n_pw}" == S20.loc[("pairwise", "tcr70: TCR-IMAP vs pair-frequency prior"), "wins"]),
        ("win rate vs TEIM-Res == step-1", f"{win_teim}/{n_pw}" == S20.loc[("pairwise", "tcr70: TCR-IMAP vs TEIM-Res"), "wins"]),
        ("position prior (tcr70, complex-macro) CI includes zero",
         not ci_excludes_zero(perres[("tcr70", "PosFreq")]["lo"], perres[("tcr70", "PosFreq")]["hi"]))]
print("\nstep-1 (stats_20k_seed42.tsv) cross-check:")
for name, ok in chk:
    print(f"  {'OK      ' if ok else 'MISMATCH'} {name}")
assert all(ok for _, ok in chk)
print(f"\n  note: cluster-macro 20k p: vs Struct-TCS {ptab(macro['Struct-TCS']['p'], N_BOOT)}, "
      f"vs position prior {ptab(macro['PosFreq']['p'], N_BOOT)}, vs pair prior {ptab(macro['pair']['p'], N_BOOT)}")
print("  note: panel F 6 A single random split (20k): p = %s, CI [%+.4f, %+.4f]" % (fvals[2][6], fvals[2][4], fvals[2][5]))
n3, arial, bases = pdf_font_check(pdf)
print(f"\nPDF font check {pdf}: Type3 = {n3}, ArialMT present = {arial}, BaseFonts = {bases}")
assert n3 == 0 and arial
