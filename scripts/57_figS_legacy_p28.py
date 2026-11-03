#!/usr/bin/env python3
"""SI figure S7 - the legacy 28-position results beside the reported p84v2 values, with the selection bias disclosed.

The 28 MHC positions of the legacy analysis were chosen from contact frequency pooled over all 174 complexes BEFORE the
data were split - outcome-informed selection (HANDOFF_1.md section 2; CORRECTION_PLAN_P84V2.md section 1). These
results are superseded and shown only to document what the corrections changed: first the structurally defined
82-position window (p82v1), then the residue representation (p84v2), and, for the per-residue models, the
cluster-disjoint validation split.

Five headline comparisons, TCR-IMAP minus the comparator, under tcr70 unless stated: pairwise vs the pair-frequency
prior; per-residue vs Struct-TCS; per-residue vs the position prior; the identity sweep at < 0.70 vs the pair-frequency
prior; peptide input, CDR3 rows. Difference with 95 % CI and p, and the two macro-AUPRC values.

Sources (nothing is trained, predicted or re-scored):
- p28, as recorded in the committed August tables (2,000 resamples, seed 0), all written from
  data/processed/pair_maps.npz (174 x 84 x 28, checked here): protocol_main.csv (24_run_all.py),
  benchmark_v2.csv (28_benchmark_v2.py), similarity_curve.csv (25_similarity_curve.py), peptide_aware.csv
  (34_peptide_aware.py)
- p84v2, as reported: results/review/figures_p84v2/stats_20k_seed42.tsv for the three rows with stored per-complex
  arrays (recomputed here at 20,000 resamples, seed 42, and asserted equal); similarity_curve_p84v2.csv and
  peptide_aware_p84v2.csv (recorded 2,000 resamples, seed 0; no per-complex arrays) for the other two
Fill = 95 % CI excludes zero. Marker = analysis: diamond p28 (superseded), circle p84v2 (reported).
"""
import sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, "src")
from viz.style import (apply_style, save_fig, print_validation_table, pstr, sig_marker, signed, minus,
                       PAGE_WIDTH_IN, MODEL_LABEL, METHOD_COLORS)
from viz.figstats import pooled_cluster_bootstrap, ci_excludes_zero, parse_cell, pdf_font_check, N_BOOT, SEED

apply_style()
T = "results/tables/"
OUT_STEM = "figures/figS_legacy_p28_vs_p84v2_v1"
TXT = "0.30"
RD = dict(float_precision="round_trip")
VAL = []
def val(what, shown, src, cell):
    VAL.append(dict(what=what, shown=shown, src=src, cell=cell))

L = np.load("data/processed/pair_maps.npz", allow_pickle=True)
assert L["Y"].shape == (174, 84, 28), L["Y"].shape          # the 28-position candidate set
for f, s in (("protocol_main.csv", "24_run_all.py"), ("benchmark_v2.csv", "28_benchmark_v2.py"),
             ("similarity_curve.csv", "25_similarity_curve.py"), ("peptide_aware.csv", "34_peptide_aware.py")):
    src = open(f"scripts/{s}").read()
    assert f in src and "pair_maps.npz" in src and "pair_maps_p8" not in src, (f, s)   # written from the p28 labels
print("p28 labels are 174 x 84 x 28; each p28 table is written by a script that reads data/processed/pair_maps.npz")

# ------------------------------------------------------------------ p28, recorded
pm = pd.read_csv(T + "protocol_main.csv", **RD).set_index("protocol").loc["tcr70"]
bv = pd.read_csv(T + "benchmark_v2.csv", **RD).set_index("protocol").loc["tcr70"]
sc = pd.read_csv(T + "similarity_curve.csv", **RD).set_index("threshold").loc[0.7]
pa = pd.read_csv(T + "peptide_aware.csv", **RD).set_index("scope").loc["CDR3"]
P28 = {
    "pair":   dict(diff=pm["diff"], lo=pm.ci_lo, hi=pm.ci_hi, p=pm.pval, a=pm.model_auprc, b=pm.freq_auprc,
                   src="protocol_main.csv", cell="protocol=tcr70: diff, ci_lo, ci_hi, pval, model_auprc, freq_auprc"),
    "struct": dict(zip(("diff", "lo", "hi", "p"), parse_cell(bv["vs_Struct-TCS"])), a=bv.Ours, b=bv["Struct-TCS"],
                   src="benchmark_v2.csv", cell="protocol=tcr70: vs_Struct-TCS, Ours, Struct-TCS"),
    "posfreq": dict(zip(("diff", "lo", "hi", "p"), parse_cell(bv["vs_PosFreq"])), a=bv.Ours, b=bv.PosFreq,
                    src="benchmark_v2.csv", cell="protocol=tcr70: vs_PosFreq, Ours, PosFreq"),
    "sweep":  dict(diff=sc["diff"], lo=sc.lo, hi=sc.hi, p=sc.p, a=sc.model, b=sc.freq,
                   src="similarity_curve.csv", cell="threshold=0.7: diff, lo, hi, p, model, freq"),
    "pep":    dict(diff=pa["diff"], lo=pa.lo, hi=pa.hi, p=pa.pval, a=pa.pep, b=pa.no_pep,
                   src="peptide_aware.csv", cell="scope=CDR3: diff, lo, hi, pval, pep, no_pep"),
}

# ------------------------------------------------------------------ p84v2, as reported
S20 = pd.read_csv("results/review/figures_p84v2/stats_20k_seed42.tsv", sep="\t", **RD)
z = np.load(T + "percomplex_tcr70_p84v2.npz", allow_pickle=True)
pr = pd.read_csv(T + "percomplex_tcr70_paired_p84v2.tsv", sep="\t", **RD)
b84 = pd.read_csv(T + "benchmark_p84v2.csv", **RD).set_index("protocol").loc["tcr70"]
P84 = {}
for key, (a, b, cl), grp, name, ab in (
        ("pair", (pr.grip.values, pr.prior.values, pr.clu70.values), "pairwise", "tcr70: TCR-IMAP vs pair-frequency prior",
         (0.3598, 0.3090)),
        ("struct", (z["Ours_auprc"], z["Struct-TCS_auprc"], z["clu70"]), "per-residue", "tcr70: TCR-IMAP vs Struct-TCS",
         (b84.Ours, b84["Struct-TCS"])),
        ("posfreq", (z["Ours_auprc"], z["PosFreq_auprc"], z["clu70"]), "per-residue", "tcr70: TCR-IMAP vs PosFreq",
         (b84.Ours, b84.PosFreq))):
    r = pooled_cluster_bootstrap(a, b, cl, N_BOOT, SEED)
    rec = S20[(S20.group == grp) & (S20.comparison == name)].iloc[0]
    for f in ("diff", "lo", "hi", "p"):
        assert abs(r[f] - rec[f]) < 1e-12, (key, f)
    assert round(float(np.nanmean(a)), 4) == round(ab[0], 4) and round(float(np.nanmean(b)), 4) == round(ab[1], 4), key
    P84[key] = dict(diff=r["diff"], lo=r["lo"], hi=r["hi"], p=r["p"], a=float(np.nanmean(a)), b=float(np.nanmean(b)),
                    src="stats_20k_seed42.tsv", cell=f"{grp}: {name}", boot=f"{N_BOOT // 1000}k")
sc84 = pd.read_csv(T + "similarity_curve_p84v2.csv", **RD).set_index("threshold").loc[0.7]
pa84 = pd.read_csv(T + "peptide_aware_p84v2.csv", **RD).set_index("scope").loc["CDR3"]
P84["sweep"] = dict(diff=sc84["diff"], lo=sc84.lo, hi=sc84.hi, p=sc84.p, a=sc84.model, b=sc84.freq,
                    src="similarity_curve_p84v2.csv", cell="threshold=0.7", boot="2k")
P84["pep"] = dict(diff=pa84["diff"], lo=pa84.lo, hi=pa84.hi, p=pa84.pval, a=pa84.pep, b=pa84.no_pep,
                  src="peptide_aware_p84v2.csv", cell="scope=CDR3", boot="2k")
# the CLAUDE.md headline values, as a guard against reading the wrong row
assert (round(P84["pair"]["diff"], 4), round(P84["struct"]["diff"], 4), round(P84["posfreq"]["diff"], 4),
        round(P84["sweep"]["diff"], 4), round(P84["pep"]["diff"], 4)) == (0.0507, 0.0312, 0.0085, 0.0274, 0.0082)
print(f"p84v2: 3 rows recomputed at {N_BOOT:,} resamples, seed {SEED} = stats_20k_seed42.tsv; sweep and peptide rows "
      "recorded (2,000 resamples)")

ROWS = [("pair", "pairwise, tcr70\nvs pair-frequency prior", METHOD_COLORS["Pair-frequency prior"], "pair prior"),
        ("struct", "per-residue, tcr70\nvs Struct-TCS", METHOD_COLORS["Struct-TCS"], "Struct-TCS"),
        ("posfreq", "per-residue, tcr70\nvs position prior", METHOD_COLORS["Position prior"], "position prior"),
        ("sweep", "identity sweep, < 0.70\nvs pair-frequency prior", METHOD_COLORS["Pair-frequency prior"], "pair prior"),
        ("pep", "peptide input, CDR3 rows\nwith vs without", METHOD_COLORS[MODEL_LABEL], "without")]
VER = {"p28": ("D", "p28: 28 positions, outcome-informed (superseded)"), "p84v2": ("o", "p84v2: 84-column window (reported)")}

# ------------------------------------------------------------------ figure
FW, FH = PAGE_WIDTH_IN, 3.55
fig = plt.figure(figsize=(FW, FH))
LA, WA = 1.62, 1.70
ax = fig.add_axes([LA / FW, 0.88 / FH, WA / FW, 2.30 / FH])
XD, XP, XM, XB = 1 + 0.10 / WA, 1 + 1.42 / WA, 1 + 2.05 / WA, 1 + 3.23 / WA
yy = np.arange(len(ROWS))[::-1].astype(float)
for (key, lab, col, blab), y in zip(ROWS, yy):
    for dy, v, R in ((+0.17, "p28", P28[key]), (-0.17, "p84v2", P84[key])):
        sig = ci_excludes_zero(R["lo"], R["hi"])
        ax.plot([R["lo"], R["hi"]], [y + dy, y + dy], color=col, lw=1.0, solid_capstyle="butt")
        ax.plot(R["diff"], y + dy, VER[v][0], ms=4.2 if v == "p28" else 4.6, zorder=3, **sig_marker(col, sig))
        tc = "0.55" if v == "p28" else TXT
        ax.text(XD, y + dy, minus(f"{signed(R['diff'], 4)} [{signed(R['lo'], 4)}, {signed(R['hi'], 4)}]"),
                transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.7, color=tc)
        ax.text(XP, y + dy, pstr(R["p"]), transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.7, color=tc)
        ax.text(XM, y + dy, f"{R['a']:.4f} vs {R['b']:.4f}", transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=5.7, color=tc)
        boot = "2k" if v == "p28" else R["boot"]
        ax.text(XB, y + dy, f"{v}, {boot}", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.7, color=tc)
        val(f"{lab.replace(chr(10), ' ')} [{v}]", f"{signed(R['diff'], 4)} [{signed(R['lo'], 4)}, {signed(R['hi'], 4)}] "
            f"{pstr(R['p'])}; {R['a']:.4f} vs {R['b']:.4f} ({'filled' if sig else 'hollow'}; {boot})", R["src"], R["cell"])
ax.axvline(0, ls="--", lw=0.8, color="0.4")
ax.set_yticks(yy); ax.set_yticklabels([r[1] for r in ROWS], fontsize=6.2, linespacing=1.15)
ax.set_ylim(-0.6, len(ROWS) - 0.4)
ax.set_xlabel(f"{MODEL_LABEL} − comparator (macro-AUPRC)", fontsize=7)
for x, h in ((XD, "Δ [95% CI]"), (XP, "p"), (XM, f"{MODEL_LABEL} vs comparator"), (XB, "resamples")):
    ax.text(x, 1.0, h, transform=ax.transAxes, ha="left", va="bottom", fontsize=6.2, color="0.2")
fig.legend(handles=[Line2D([], [], marker=m, ls="", mfc="0.3", mec="0.3", ms=4.3, label=l) for m, l in VER.values()],
           loc="lower left", bbox_to_anchor=(0.01, 0.075), ncol=2, frameon=False, fontsize=6.3, handletextpad=0.3,
           columnspacing=1.5)
ax.set_title("Legacy 28-position results beside p84v2", fontsize=8, pad=12, loc="left", x=-0.85)
fig.text(0.01, 0.012, minus("p28 positions were selected from contact frequency pooled over all 174 complexes before splitting "
                            "(outcome-informed); superseded, shown only to document the correction"),
         ha="left", va="bottom", fontsize=5.6, color="0.25")
png, pdf = save_fig(fig, OUT_STEM)
plt.close(fig)

print("\n--- every drawn number and its source ---")
print_validation_table(VAL, [("what", "row", None), ("shown", "drawn", None), ("src", "source", None), ("cell", "cell", None)])
n3, arial, bases = pdf_font_check(pdf)
print(f"\nPDF font check {pdf}: Type3 = {n3}, ArialMT present = {arial}, BaseFonts = {bases}")
assert n3 == 0 and arial
