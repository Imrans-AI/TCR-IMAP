#!/usr/bin/env python3
"""SI figure (p84v2) - peptide information: no detectable effect.

Built from 34c_figS_peptide_p82v1.py; input results/tables/peptide_aware_p84v2.csv
(scripts/53_peptide_aware_p84v2.py). Under p84v2 every row is null: no CI excludes zero, CDR3
+0.0082 [-0.0137, +0.0293], p = 0.416 (p82v1: +0.0234, p = 0.053). Title and legend say so.

TCR-IMAP with vs without peptide input, tcr70, pairwise macro-AUPRC difference with 95 %
cluster-bootstrap CI: overall, shuffled-peptide control, germline loops (CDR1/2 rows only) and
CDR3 rows only (rows split by the p84v2 row segments).

These are the RECORDED values (2,000 resamples, seed 0). No per-complex AUPRC arrays exist on
disk for this analysis -- the 34b cache and peptide_preds_p84v2.npz hold predictions only, and
scoring them would be re-evaluation, which this figure pass does not do. So no 20k recompute.

(scripts/55_stats_p84v2.py did score the stored predictions for provenance: they reproduce every
recorded value exactly at 2,000 resamples, and its 20,000-resample values are in
results/review/figures_p84v2/stats_20k_seed42.tsv; they are not displayed here.)
Framing: no peptide effect is detected in any row; SI only, not evidence for two-step docking.
Fill = 95 % CI excludes zero (none does here); p via style.pstr() only.
"""
import sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, "src")
from viz.style import (apply_style, save_fig, print_validation_table, pstr, sig_marker, signed, minus,
                       PAGE_WIDTH_IN, MODEL_LABEL, METHOD_COLORS)
from viz.figstats import ci_excludes_zero, pdf_font_check

apply_style()
OUT_STEM = "figures/figS_peptide_p84v2_v1"
CM = METHOD_COLORS[MODEL_LABEL]
TXT = "0.30"

pep = pd.read_csv("results/tables/peptide_aware_p84v2.csv").set_index("scope")
ROWS = [("peptide", "all positions"),
        ("peptide_shuf", "all positions, shuffled-peptide control"),
        ("germline(CDR1/2)", "CDR1/2 rows only"),
        ("CDR3", "CDR3 rows only")]
assert list(pep.index) == [r[0] for r in ROWS]

fig = plt.figure(figsize=(PAGE_WIDTH_IN, 2.35))
FW, FH = PAGE_WIDTH_IN, 2.35
ax = fig.add_axes([2.05 / FW, 0.48 / FH, 2.00 / FW, 1.45 / FH])
yy = np.arange(len(ROWS))[::-1]
for (key, lab), y in zip(ROWS, yy):
    r = pep.loc[key]
    ax.plot([r.lo, r.hi], [y, y], color=CM, lw=1.1, solid_capstyle="butt")
    ax.plot(r["diff"], y, "o", ms=4.6, zorder=3, **sig_marker(CM, ci_excludes_zero(r.lo, r.hi)))
    ax.text(1.03, y, f"{signed(r['diff'], 4)} [{signed(r.lo, 4)}, {signed(r.hi, 4)}]", transform=ax.get_yaxis_transform(),
            ha="left", va="center", fontsize=6, color=TXT)
    ax.text(1.78, y, pstr(r.pval), transform=ax.get_yaxis_transform(), ha="left", va="center",
            fontsize=6, color=TXT)
    ax.text(2.12, y, minus(f"{r.no_pep:.4f} → {r.pep:.4f}"), transform=ax.get_yaxis_transform(),
            ha="left", va="center", fontsize=6, color=TXT)
for x, h in ((1.03, "Δ [95% CI]"), (1.78, "p"), (2.12, "without → with")):
    ax.text(x, 1.02, h, transform=ax.transAxes, ha="left", va="bottom", fontsize=6.2, color="0.2")
ax.axvline(0, ls="--", lw=0.8, color="0.4")
ax.set_yticks(yy); ax.set_yticklabels([lab for _, lab in ROWS], fontsize=6.5)
ax.set_ylim(-0.6, len(ROWS) - 0.4)
ax.set_xlim(min(-0.02, pep.lo.min() - 0.005), max(0.04, pep.hi.max() + 0.008))
ax.set_xlabel(f"{MODEL_LABEL} with − without peptide (pairwise macro-AUPRC)", fontsize=7)
assert not any(ci_excludes_zero(pep.loc[k].lo, pep.loc[k].hi) for k, _ in ROWS)
ax.set_title("Peptide input, tcr70: no detectable effect (every 95% CI includes zero)", fontsize=8, pad=12, loc="left",
             x=-0.92)
fig.text(0.99, 0.02, minus("recorded values, clu70 cluster bootstrap, 2,000 resamples"), ha="right",
         va="bottom", fontsize=5.4, color="0.45")

png, pdf = save_fig(fig, OUT_STEM)
plt.close(fig)

rows = [dict(scope=k, no_pep=f"{pep.loc[k].no_pep:.4f}", pep=f"{pep.loc[k].pep:.4f}",
             diff=f"{pep.loc[k]['diff']:+.4f}", ci=f"[{pep.loc[k].lo:+.4f}, {pep.loc[k].hi:+.4f}]",
             p=f"{pep.loc[k].pval:.3f}", shown=pstr(pep.loc[k].pval),
             fill="filled" if ci_excludes_zero(pep.loc[k].lo, pep.loc[k].hi) else "hollow")
        for k, _ in ROWS]
print("\n--- peptide_aware_p84v2.csv (recorded, 2k; no per-complex arrays -> no 20k recompute) ---")
print_validation_table(rows, [("scope", "scope", None), ("no_pep", "without", None), ("pep", "with", None),
                              ("diff", "diff", None), ("ci", "95% CI", None), ("p", "p recorded", None),
                              ("shown", "displayed", None), ("fill", "marker", None)])
EXP = {"peptide": (0.0100, 0.196), "peptide_shuf": (0.0071, 0.308),
       "germline(CDR1/2)": (0.0107, 0.133), "CDR3": (0.0082, 0.416)}
print("\ncross-check against the step-1 provenance table (recorded p84v2, 2k):")
for k, (d, p) in EXP.items():
    r = pep.loc[k]
    ok = abs(round(r["diff"], 4) - d) < 1e-9 and abs(r.pval - p) < 1e-9
    print(f"  {'OK      ' if ok else 'MISMATCH'} {k}: {r['diff']:+.4f} p {r.pval:.3f} (expected {d:+.4f} / {p})")
    assert ok
S20 = pd.read_csv("results/review/figures_p84v2/stats_20k_seed42.tsv", sep="\t")
print("\nfor reference only (not displayed): 20,000-resample values from stored predictions (step 1):")
for _, r in S20[S20.group == "peptide"].iterrows():
    print(f"  {r.comparison:42s} {r['diff']:+.4f} [{r.lo:+.4f}, {r.hi:+.4f}] p {r.p:.4f}")
n3, arial, bases = pdf_font_check(pdf)
print(f"\nPDF font check {pdf}: Type3 = {n3}, ArialMT present = {arial}, BaseFonts = {bases}")
assert n3 == 0 and arial
