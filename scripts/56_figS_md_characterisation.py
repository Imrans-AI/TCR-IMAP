#!/usr/bin/env python3
"""SI figure - MD simulation stability and interface characterisation, per system and replicate. Descriptive only.

Reads the characterisation tables as recomputed with periodic-boundary handling (results/tables/
md_characterisation_*.csv, md_rmsf_per_residue.csv, md_hbond_pairs.csv; md_characterisation_manifest.csv; criteria
document, "Characterisation recomputed with periodic-boundary handling", 24 October 2026). Nothing is recomputed.

Rows (columns = the four systems, replicates coloured by REPLICATE_COLORS):
A  interface RMSD (CDR-loop Ca after a fit on MHC alpha1/alpha2 Ca, vs the crystal pose), every frame
B  buried SASA of the TCR-pMHC interface, every 1 ns
C  per-loop Ca RMSF of the canonical CDR loops (summary table), loop labels coloured by TCR chain
D  interface hydrogen bonds per frame, TCR-MHC alpha1/alpha2 and TCR-peptide: grouped bars (mean) with SD whiskers,
   replicates side by side - no per-frame H-bond series exists in any committed table, so no time course is possible
Backbone RMSD moved to Figure 7 panel E (24 October 2026). The tables hold no radius of gyration.

Cross-checks printed and asserted: the manifest confirms the PBC-handled recompute and that the occupancy tables are
unchanged; every summary mean drawn here (interface RMSD, buried SASA) equals the mean of its time series; every per-loop RMSF in
the summary equals the mean of its loop's rows in md_rmsf_per_residue.csv.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

sys.path.insert(0, "."); sys.path.insert(0, "src")
from viz.style import (apply_style, save_fig, print_validation_table, PAGE_WIDTH_IN, PANEL_LABEL_STYLE as PL,
                       REPLICATE_COLORS, MOLECULE_COLORS, panel_letters)
from viz.figstats import pdf_font_check

apply_style()
T = "results/tables/"
OUT_STEM = "figures/figS_md_characterisation_v1"
SYSTEMS = ["5yxn", "7n1f", "3vxs", "8eo8"]
REPS = (1, 2, 3)
LOOPS = ["cdr1a", "cdr2a", "cdr3a", "cdr1b", "cdr2b", "cdr3b"]
LOOP_LAB = {"cdr1a": "1α", "cdr2a": "2α", "cdr3a": "3α", "cdr1b": "1β", "cdr2b": "2β", "cdr3b": "3β"}
TXT = "0.30"
RD = dict(float_precision="round_trip")

man = pd.read_csv(T + "md_characterisation_manifest.csv").set_index("item")["value"]
assert man["occupancy_tables_unchanged"] == "True" and "nojump" in man["coordinate_measures_trajectory"]
summ = pd.read_csv(T + "md_characterisation_summary.csv", **RD)
ts = pd.read_csv(T + "md_characterisation_timeseries.csv", **RD)
rmsf = pd.read_csv(T + "md_rmsf_per_residue.csv", **RD)
print(f"characterisation recomputed at {man['git_commit'][:7]} with periodic-boundary handling; occupancy tables unchanged")

# ------------------------------------------------------------------ cross-checks
XC = []
for p in SYSTEMS:
    for r in REPS:
        s = summ[(summ.pdb == p) & (summ.rep == r)].iloc[0]
        t = ts[(ts.pdb == p) & (ts.rep == r)]
        for col, scol in (("interface_rmsd", "interface_rmsd_mean"), ("buried_sasa", "buried_sasa_mean")):
            v = t[col].dropna()
            assert np.isclose(v.mean(), s[scol], rtol=0, atol=1e-9), (p, r, col)
        for L in LOOPS:
            m = rmsf[(rmsf.pdb == p) & (rmsf.rep == r) & (rmsf.loop == L)].rmsf.mean()
            assert np.isclose(m, s[f"rmsf_{L}"], rtol=0, atol=1e-9), (p, r, L)
        XC.append(dict(sys=f"{p} rep{r}", frames=int(t.interface_rmsd.notna().sum()), sasa=int(t.buried_sasa.notna().sum()),
                       ifr=f"{s.interface_rmsd_mean:.2f}", sasa_m=f"{s.buried_sasa_mean:.0f}",
                       hb=f"{s['hbond_TCR-MHC_mean']:.2f} / {s['hbond_TCR-peptide_mean']:.2f}"))
print("\nsummary means = means of their time series; per-loop RMSF = mean of per-residue rows (all asserted):")
print_validation_table(XC, [("sys", "replicate", None), ("frames", "RMSD frames", None), ("sasa", "SASA frames", None),
                            ("ifr", "interface RMSD mean (A)", None),
                            ("sasa_m", "buried SASA mean (A^2)", None), ("hb", "H-bonds TCR-MHC / TCR-pep", None)])

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(PAGE_WIDTH_IN, 7.4))
gs = fig.add_gridspec(4, 4, hspace=0.68, wspace=0.12, left=0.085, right=0.975, top=0.95, bottom=0.085,
                      height_ratios=[1, 1, 0.95, 0.85])
rows = {}
YL = {"interface_rmsd": (0, 12.0), "buried_sasa": (0, 2400)}
YLAB = {"interface_rmsd": "Interface RMSD (Å)", "buried_sasa": "Buried SASA (Å²)"}
for ri, col in enumerate(["interface_rmsd", "buried_sasa"]):
    rows[ri] = []
    for ci, p in enumerate(SYSTEMS):
        ax = fig.add_subplot(gs[ri, ci]); rows[ri].append(ax)
        for r in REPS:
            t = ts[(ts.pdb == p) & (ts.rep == r)].dropna(subset=[col])
            ax.plot(t.time_ps / 1000, t[col], color=REPLICATE_COLORS[r], lw=0.45 if col != "buried_sasa" else 0.7,
                    rasterized=True, zorder=4 - r)
        if col == "buried_sasa":                # thousands separators on the SASA axis
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.set_ylim(*YL[col]); ax.set_xlim(20, 200); ax.set_xticks([50, 100, 150, 200])   # analysed window, to 200 ns
        if ci == 0:
            ax.set_ylabel(YLAB[col])
        else:
            ax.tick_params(labelleft=False)
        if ri == 0:
            ax.set_title(p, fontsize=8, pad=3)
        ax.set_xlabel("Time (ns)")                 # every time-axis row carries its own axis

# C: per-loop RMSF
rows[2] = []
for ci, p in enumerate(SYSTEMS):
    ax = fig.add_subplot(gs[2, ci]); rows[2].append(ax)
    for li, L in enumerate(LOOPS):
        for r in REPS:
            v = summ[(summ.pdb == p) & (summ.rep == r)][f"rmsf_{L}"].iloc[0]
            x = li + (r - 2) * 0.22
            ax.plot([x - 0.09, x + 0.09], [v, v], color=REPLICATE_COLORS[r], lw=1.6, solid_capstyle="butt")
    ax.set_xticks(range(len(LOOPS))); ax.set_xticklabels([LOOP_LAB[L] for L in LOOPS], fontsize=6.5)
    for lab, L in zip(ax.get_xticklabels(), LOOPS):
        lab.set_color(MOLECULE_COLORS["TCR alpha" if L.endswith("a") else "TCR beta"])
    ax.set_ylim(0, 4.5); ax.set_xlim(-0.6, len(LOOPS) - 0.4)
    if ci == 0:
        ax.set_ylabel("Loop Cα RMSF (Å)")
    else:
        ax.tick_params(labelleft=False)
    ax.set_xlabel("CDR loop", fontsize=7)

# D: H-bonds
rows[3] = []
for ci, p in enumerate(SYSTEMS):
    ax = fig.add_subplot(gs[3, ci]); rows[3].append(ax)
    # No committed table holds a per-frame H-bond series (only the window mean and SD per replicate, and per-residue-
    # pair occupancy), so this is drawn as grouped bars: mean per frame, whisker = frame-to-frame SD over 20-200 ns,
    # lower whisker clipped at 0 H-bonds. Replicates side by side; bar colour = replicate, nothing else.
    for ii, tag in enumerate(["TCR-MHC", "TCR-peptide"]):
        for r in REPS:
            s = summ[(summ.pdb == p) & (summ.rep == r)].iloc[0]
            m, sd = s[f"hbond_{tag}_mean"], s[f"hbond_{tag}_sd"]
            x = ii + (r - 2) * 0.27
            ax.bar(x, m, 0.25, color=REPLICATE_COLORS[r], edgecolor="none", zorder=2)
            ax.errorbar(x, m, yerr=[[min(sd, m)], [sd]], fmt="none", ecolor="0.15", elinewidth=0.6, capsize=1.6,
                        capthick=0.6, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["TCR–MHC", "TCR–peptide"], fontsize=6.5)
    ax.set_ylim(0, 8.5); ax.set_xlim(-0.5, 1.5)
    if ci == 0:
        ax.set_ylabel("H-bonds per frame")
    else:
        ax.tick_params(labelleft=False)

fig.legend(handles=[Line2D([], [], color=REPLICATE_COLORS[r], lw=1.6, label=f"replicate {r}") for r in REPS],
           loc="lower center", ncol=3, frameon=False, fontsize=6.5, bbox_to_anchor=(0.54, 0.0))

panel_letters(fig, {letter: rows[ri] for ri, letter in enumerate("ABCD")})   # one fixed offset (style.py)
png, pdf = save_fig(fig, OUT_STEM)
n3, arial, bases = pdf_font_check(pdf)
print(f"\nPDF font check {pdf}: Type3 = {n3}, ArialMT present = {arial}, BaseFonts = {bases}")
assert n3 == 0 and arial
