#!/usr/bin/env python3
"""SI figure S4 (p84v2) - sensitivity of the per-residue benchmark to the validation split used for model selection.

The per-residue models (TCR-IMAP's per-residue variant, Struct-TCS, Seq-TCS) early-stop on a validation subset of the
training complexes. The reported analysis draws it as whole clu70 clusters (cluster-disjoint); the superseded analysis
drew a random 15 % of complexes. Predeclared in CORRECTION_PLAN_P84V2.md; full record VALIDATION_FIX_REPORT.md.

A  per-residue macro-AUPRC under both validation splits, every protocol, the three trained models; the change in
   TCR-IMAP's macro-AUPRC is written above each protocol
B  tcr70: TCR-IMAP vs Struct-TCS, vs Seq-TCS and vs the position prior under both splits, 95 % CI and p
C  clu70 clusters with complexes on both sides of model selection, per fold, under the random split (18 folds);
   zero in every fold under the cluster-disjoint split

Sources (nothing is trained, predicted or re-scored):
- results/review/validation_fix/randomsplit_deterministic/percomplex_{protocol}_p84v2.npz (random split; sha256
  checked against its committed SHA256SUMS) and results/tables/percomplex_{protocol}_p84v2.npz (cluster-disjoint;
  checked against clusterdisjoint_run1/SHA256SUMS)
- results/review/validation_fix/B3_macro_auprc.tsv (macro-AUPRC, both splits) - every drawn mean is asserted equal to
  the mean of the stored per-complex array
- results/review/validation_fix/B3_key_differences_20k.tsv (paired differences, pooled clu70 cluster bootstrap, 20,000
  resamples, seed 42) - recomputed here from the stored arrays with figstats.pooled_cluster_bootstrap and asserted equal
- results/review/validation_fix/queue_console_randomsplit_deterministic.log (committed run log): the per-fold count of
  clu70 clusters on both sides of model selection, printed by the sensitivity run itself. No table holds these counts;
  they are read from that log, not recomputed. results/review/validation_fix/B2_split_dryrun.log records that every
  cluster-disjoint fold is disjoint.
Fill = 95 % CI excludes zero (panel B only). Marker shape = validation split: triangle random, circle cluster-disjoint.
"""
import sys, re, hashlib
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, "src")
from viz.style import (apply_style, save_fig, print_validation_table, pstr, sig_marker, signed, minus, panel_letters,
                       PAGE_WIDTH_IN, MODEL_LABEL, METHOD_COLORS, PROTOCOL_COLORS)
from viz.figstats import pooled_cluster_bootstrap, ci_excludes_zero, pdf_font_check, N_BOOT, SEED

apply_style()
T, V = "results/tables/", "results/review/validation_fix/"
OUT_STEM = "figures/figS_validation_split_p84v2_v1"
PROTOS = ["random", "double", "allele", "temporal", "tcr70"]       # light -> dark strictness
TRAINED = {"Ours": MODEL_LABEL, "Struct-TCS": "Struct-TCS", "Seq-TCS": "Seq-TCS"}
COMP = {"Struct-TCS": "Struct-TCS", "Seq-TCS": "Seq-TCS", "PosFreq": "Position prior"}
SPLITS = {"random split": V + "randomsplit_deterministic/", "cluster-disjoint": T}
SHA = {"random split": V + "randomsplit_deterministic/SHA256SUMS", "cluster-disjoint": V + "clusterdisjoint_run1/SHA256SUMS"}
MK = {"random split": "^", "cluster-disjoint": "o"}
TXT = "0.30"
RD = dict(float_precision="round_trip")
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
VAL = []
def val(panel, what, shown, src, cell):
    VAL.append(dict(panel=panel, what=what, shown=shown, src=src, cell=cell))

# ------------------------------------------------------------------ inputs and cross-checks
pcx = {}
for s, d in SPLITS.items():
    sums = dict(l.split()[::-1] for l in open(SHA[s]) if l.strip())
    for p in PROTOS:
        f = f"percomplex_{p}_p84v2.npz"
        assert sha(d + f) == sums[f], (s, f)                # exactly the arrays the committed checksums describe
        pcx[(s, p)] = np.load(d + f, allow_pickle=True)
    assert sha(d + "benchmark_p84v2.csv") == sums["benchmark_p84v2.csv"], s
print("all 10 per-complex arrays and both benchmark tables match their committed SHA256SUMS")

mac = pd.read_csv(V + "B3_macro_auprc.tsv", sep="\t", **RD)
KEY = {"Ours": "TCR-IMAP", "Struct-TCS": "Struct-TCS", "Seq-TCS": "Seq-TCS", "PosFreq": "Position prior", "NN": "NN transfer"}
M = {}
for p in PROTOS:
    for k, lab in KEY.items():
        r = mac[(mac.protocol == p) & (mac.method == lab)].iloc[0]
        for s, col in (("random split", "random_split"), ("cluster-disjoint", "cluster_disjoint")):
            z = pcx[(s, p)]
            assert (z["pdb"] == pcx[("cluster-disjoint", p)]["pdb"]).all()     # same test complexes, same order
            m = float(np.nanmean(z[f"{k}_auprc"]))
            assert abs(m - r[col]) < 1e-12, (p, k, s, m, r[col])
            M[(s, p, k)] = m
        if k in ("PosFreq", "NN"):                         # fitted without a validation subset: must not move at all
            assert np.array_equal(pcx[("random split", p)][f"{k}_auprc"], pcx[("cluster-disjoint", p)][f"{k}_auprc"],
                                  equal_nan=True), (p, k)
print("macro-AUPRC: every cell of B3_macro_auprc.tsv equals the mean of its stored array; position prior and NN "
      "transfer arrays identical under both splits (5/5 protocols)")

kd = pd.read_csv(V + "B3_key_differences_20k.tsv", sep="\t", **RD)
D = {}
for s in SPLITS:
    z = pcx[(s, "tcr70")]
    for k, lab in COMP.items():
        r = pooled_cluster_bootstrap(z["Ours_auprc"], z[f"{k}_auprc"], z["clu70"], N_BOOT, SEED)
        rec = kd[(kd.protocol == "tcr70") & (kd.comparison == f"TCR-IMAP vs {lab}") & (kd.split == s)].iloc[0]
        for a, b in (("diff", "diff"), ("lo", "lo"), ("hi", "hi"), ("p", "p")):
            assert abs(r[a] - rec[b]) < 1e-12, (s, k, a, r[a], rec[b])
        D[(s, k)] = r
print(f"tcr70 key differences recomputed at {N_BOOT:,} resamples, seed {SEED}: 6/6 equal B3_key_differences_20k.tsv")

# straddling clusters per fold, from the sensitivity run's own log
STR = []
pat = re.compile(r"\[(\w+)\] SENSITIVITY random split: test (\d+), train (\d+), validation (\d+) complexes = ([\d.]+)% "
                 r"of (\d+) training complexes, NOT cluster-disjoint \((\d+) clu70 clusters on both sides\)")
for l in open(V + "queue_console_randomsplit_deterministic.log"):
    m = pat.search(l)
    if m:
        STR.append(dict(protocol=m[1], test=int(m[2]), val=int(m[4]), straddle=int(m[7])))
STR = pd.DataFrame(STR)
STR["fold"] = STR.groupby("protocol").cumcount()
assert len(STR) == 18 and (STR.protocol.value_counts().reindex(PROTOS).values == [5, 3, 4, 1, 5]).all()
dry = open(V + "B2_split_dryrun.log").read()
assert "folds evaluated: 18; cluster-disjoint in every one: True" in dry
smin, smax, smed = int(STR.straddle.min()), int(STR.straddle.max()), float(STR.straddle.median())
print(f"straddling clusters, random split: {STR.straddle.tolist()}  min {smin}, max {smax}, median {smed}")

# ------------------------------------------------------------------ figure
FW, FH = PAGE_WIDTH_IN, 6.2
fig = plt.figure(figsize=(FW, FH))
def axin(l, b, w, h): return fig.add_axes([l / FW, b / FH, w / FW, h / FH])
C = {k: METHOD_COLORS[v] for k, v in {**TRAINED, **COMP}.items()}

# --- A
ax = axA = axin(0.62, 3.40, 6.45, 2.05)
OFF = {"Ours": -0.26, "Struct-TCS": 0.0, "Seq-TCS": 0.26}
DS = {"random split": -0.07, "cluster-disjoint": 0.07}
for i, p in enumerate(PROTOS):
    for k in TRAINED:
        xs = [i + OFF[k] + DS[s] for s in SPLITS]
        ys = [M[(s, p, k)] for s in SPLITS]
        ax.plot(xs, ys, color=C[k], lw=0.8, zorder=2)
        for s, x, y in zip(SPLITS, xs, ys):
            ax.plot(x, y, MK[s], ms=4.2, color=C[k], mec=C[k], zorder=3)
            val("A", f"{p} {TRAINED[k]} {s}", f"{y:.4f}", "B3_macro_auprc.tsv",
                f"protocol={p}, method={KEY[k]}, {'random_split' if s == 'random split' else 'cluster_disjoint'}")
    ch = M[("cluster-disjoint", p, "Ours")] - M[("random split", p, "Ours")]
    chr_ = float(mac[(mac.protocol == p) & (mac.method == "TCR-IMAP")].change.iloc[0])
    assert abs(ch - chr_) < 1e-12
    ax.text(i, 0.905, minus(f"{MODEL_LABEL} {signed(ch, 4)}"), ha="center", va="bottom", fontsize=6.2, color=TXT)
    val("A", f"{p} {MODEL_LABEL} change", signed(ch, 4), "B3_macro_auprc.tsv", f"protocol={p}, method=TCR-IMAP, change")
    if i < len(PROTOS) - 1:
        ax.axvline(i + 0.5, color="0.88", lw=0.6, zorder=0)
ax.set_xticks(range(len(PROTOS)))
ax.set_xticklabels([f"{p}\nn = {len(pcx[('cluster-disjoint', p)]['pdb'])}" + ("\n(primary)" if p == "tcr70" else "")
                    for p in PROTOS], fontsize=6.5, linespacing=1.25)
ax.set_xlim(-0.5, len(PROTOS) - 0.5); ax.set_ylim(0.60, 0.95)
ax.set_ylabel("Per-residue macro-AUPRC")
ax.text(0.5, 0.995, "change in TCR-IMAP macro-AUPRC, random split → cluster-disjoint", transform=ax.transAxes,
        ha="center", va="bottom", fontsize=6.2, color=TXT)
ax.legend(handles=[Line2D([], [], marker="s", ls="", ms=5, mfc=C[k], mec="none", label=TRAINED[k]) for k in TRAINED] +
          [Line2D([], [], marker=MK[s], ls="", ms=4.5, mfc="0.35", mec="0.35", label=f"{s} validation") for s in SPLITS],
          loc="lower left", ncol=5, frameon=False, fontsize=6.3, handletextpad=0.2, columnspacing=1.0,
          bbox_to_anchor=(0.0, -0.01))
ax.set_title("Per-residue benchmark under each validation split", fontsize=8, pad=13)

# --- B
BL, BW = 1.30, 1.45
ax = axB = axin(BL, 0.55, BW, 1.95)
XD, XP = 1 + 0.10 / BW, 1 + 1.35 / BW      # text columns at fixed inch offsets right of the axes
yy, ylab = [], []
y = 0
for k, lab in COMP.items():
    for s in SPLITS:
        r = D[(s, k)]
        ax.plot([r["lo"], r["hi"]], [y, y], color=C[k], lw=1.0, solid_capstyle="butt")
        ax.plot(r["diff"], y, MK[s], ms=4.6, zorder=3, **sig_marker(C[k], ci_excludes_zero(r["lo"], r["hi"])))
        ax.text(XD, y, minus(f"{signed(r['diff'], 4)} [{signed(r['lo'], 4)}, {signed(r['hi'], 4)}]"),
                transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.8, color=TXT)
        ax.text(XP, y, pstr(r["p"]), transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.8,
                color=TXT)
        val("B", f"tcr70 vs {lab}, {s}", f"{signed(r['diff'], 4)} [{signed(r['lo'], 4)}, {signed(r['hi'], 4)}] {pstr(r['p'])} "
            f"({'filled' if ci_excludes_zero(r['lo'], r['hi']) else 'hollow'})", "B3_key_differences_20k.tsv",
            f"protocol=tcr70, comparison=TCR-IMAP vs {lab}, split={s}")
        yy.append(y); ylab.append(f"vs {lab}\n{s}")
        y -= 1
    y -= 0.6
ax.axvline(0, ls="--", lw=0.8, color="0.4")
ax.set_yticks(yy); ax.set_yticklabels(ylab, fontsize=6.0, linespacing=1.1)
ax.set_ylim(min(yy) - 0.7, 0.7)
ax.set_xlabel(f"{MODEL_LABEL} − comparator\n(per-residue macro-AUPRC, tcr70)", fontsize=7, linespacing=1.2)
for x, h in ((XD, "Δ [95% CI]"), (XP, "p")):
    ax.text(x, 1.0, h, transform=ax.transAxes, ha="left", va="bottom", fontsize=6.2, color="0.2")
ax.set_title("tcr70 comparisons", fontsize=8, pad=12)

# --- C
ax = axC = axin(5.40, 0.55, 1.70, 1.95)
x = 0
xt, xl = [], []
for p in PROTOS:
    sub = STR[STR.protocol == p]
    xs = x + np.arange(len(sub))
    ax.bar(xs, sub.straddle, 0.75, color=PROTOCOL_COLORS[p], edgecolor="none")
    for xx, (_, r) in zip(xs, sub.iterrows()):
        val("C", f"{p} fold {r.fold} straddling clusters", str(r.straddle), "queue_console_randomsplit_deterministic.log",
            f"[{p}] SENSITIVITY random split line {r.fold + 1} of {len(sub)}")
    xt.append(xs.mean()); xl.append(p)
    x += len(sub) + 0.8
ax.axhline(smed, **dict(color="k", ls=(0, (1, 1.3)), lw=1.0))
ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=6.3, rotation=35, ha="right", rotation_mode="anchor")
ax.set_xlim(-0.8, x - 0.5); ax.set_ylim(0, 17.5)
ax.set_ylabel("clu70 clusters on both sides\nof model selection, per fold", linespacing=1.2)
ax.set_title("Random split: shared clusters", fontsize=8, pad=20)
ax.text(0.0, 1.0, minus(f"{len(STR)} folds: range {smin}–{smax}, median {smed:g} (dotted)\n"
                        "cluster-disjoint: 0 in every fold"),
        transform=ax.transAxes, ha="left", va="bottom", fontsize=5.8, color=TXT)
val("C", "range / median", f"{smin}-{smax} / {smed:g}", "queue_console_randomsplit_deterministic.log", "all 18 folds")

fig.text(0.99, 0.008, minus(f"B: pooled clu70 cluster bootstrap, {N_BOOT:,} resamples, seed {SEED}, recomputed from stored "
                            "per-complex arrays"), ha="right", va="bottom", fontsize=5.4, color="0.45")
panel_letters(fig, {"A": [axA], "B": [axB], "C": [axC]})
png, pdf = save_fig(fig, OUT_STEM)
plt.close(fig)

print("\n--- every drawn number and its source ---")
print_validation_table(VAL, [("panel", "panel", None), ("what", "value", None), ("shown", "drawn", None),
                             ("src", "source", None), ("cell", "cell", None)])
n3, arial, bases = pdf_font_check(pdf)
print(f"\nPDF font check {pdf}: Type3 = {n3}, ArialMT present = {arial}, BaseFonts = {bases}")
assert n3 == 0 and arial
