#!/usr/bin/env python3
"""Figure 2 - dataset composition and redundancy.
Sources: complex_metadata.tsv, tcr_identity_matrix.npy,
         percomplex_maxsim.tsv, leakage_audit_dual.tsv
Identity = Levenshtein.ratio over concatenated CDR1/2/3 (both chains).
Duplicate families reported under the allele-independent definition
(cdr3_pair + peptide); see Methods for the annotation-gap rationale."""
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
    "font.size":8,"axes.labelsize":8,"legend.fontsize":6.5,"xtick.labelsize":7,
    "ytick.labelsize":7,"axes.linewidth":0.6,"xtick.major.width":0.6,
    "ytick.major.width":0.6,"axes.spines.top":False,"axes.spines.right":False})
PL = dict(fontweight="bold", fontsize=10, va="bottom", ha="left")

CA, CB, RED = "#3b6ea5", "#c0504d", "#c0504d"
PC = {"random":"#c0504d","allele":"#e8a33d","double":"#7b6ca8",
      "temporal":"#4b9b8f","tcr70":"#2f5f8f"}
ORDER = ["random","allele","double","temporal","tcr70"]

m  = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t")
pc = pd.read_csv("results/tables/percomplex_maxsim.tsv", sep="\t")
au = pd.read_csv("results/tables/leakage_audit_dual.tsv", sep="\t")
S  = np.load("results/tables/tcr_identity_matrix.npy")

iu = np.triu_indices_from(S, k=1); off = S[iu]
D = S.copy(); np.fill_diagonal(D, -1); nn = D.max(axis=1)

fig = plt.figure(figsize=(7.2,5.0), constrained_layout=True)
gs  = fig.add_gridspec(2,6)

# --- A  V-gene usage --------------------------------------------------
ax = fig.add_subplot(gs[0,0:2])
ta, tb = m.trav.value_counts().head(8), m.trbv.value_counts().head(8)
labs = list(ta.index)+[""]+list(tb.index)
vals = list(ta.values)+[0]+list(tb.values)
cols = [CA]*8+["none"]+[CB]*8
y = np.arange(len(labs))[::-1]
ax.barh(y, vals, color=cols, height=0.72, linewidth=0)
ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=6.2)
ax.set_xlabel("Complexes"); ax.set_xlim(0,72)
for yy,v in zip(y,vals):
    if v: ax.text(v+1.5, yy, str(v), va="center", fontsize=5.8, color="0.3")
ax.text(44, 12.6, "TRAV", color=CA, fontsize=7.5, fontweight="bold", va="center")
ax.text(44,  3.4, "TRBV", color=CB, fontsize=7.5, fontweight="bold", va="center")
ax.set_title("V-gene usage (top 8)", fontsize=8, pad=3)
ax.text(-0.26,1.02,"A", transform=ax.transAxes, **PL)

# --- B  pairwise vs nearest-neighbour identity ------------------------
ax = fig.add_subplot(gs[0,2:4])
for v,c,lw,lab in ((off,"0.55",1.1,"All pairs (n=%s)" % f"{len(off):,}"),
                   (nn,RED,1.6,"Nearest neighbour (n=%d)" % len(nn))):
    sv = np.sort(v); ax.step(sv, np.arange(1,len(sv)+1)/len(sv), where="post",
                             color=c, lw=lw, label=lab)
ax.axvline(0.90, ls="--", lw=0.8, color="0.4")
ax.set_xlabel("TCR identity"); ax.set_ylabel("Cumulative fraction")
ax.set_xlim(0.25,1.02); ax.set_ylim(0,1.03)
ax.legend(frameon=False, fontsize=5.8, loc="upper left", handlelength=1.1)
frac_exact = float((nn >= 0.999).mean())
ax.annotate("%.0f%% have an\nexact twin" % (100*frac_exact),
            xy=(1.0, 1.0 - frac_exact/2), xytext=(0.72, 0.62),
            textcoords="axes fraction", fontsize=6, color=RED,
            ha="center", va="center", linespacing=1.4,
            arrowprops=dict(arrowstyle="-", lw=0.7, color=RED,
                            shrinkA=0, shrinkB=2))
ax.text(0.03, 0.30, "at identity \u2265 0.90:\n%.1f%% of pairs\n%.1f%% of complexes"
        % (100*(off>=0.9).mean(), 100*(nn>=0.9).mean()),
        transform=ax.transAxes, fontsize=6, va="top", ha="left",
        linespacing=1.45, color=RED)
ax.set_title("Redundancy: pairs vs neighbours", fontsize=8, pad=3)
ax.text(-0.26,1.02,"B", transform=ax.transAxes, **PL)

# --- C  cluster-size mass ---------------------------------------------
ax = fig.add_subplot(gs[0,4:6])
BL = ["1","2","3","4-5","6-10",">10"]
def b(s): return 0 if s==1 else 1 if s==2 else 2 if s==3 else 3 if s<=5 else 4 if s<=10 else 5
w = 0.26
for k,(col,cc) in enumerate(zip(["clu70","clu80","clu90"],["#2f5f8f","#6f9ac4","#b9cde0"])):
    vc = m[col].value_counts(); mass = np.zeros(6)
    for _,sz in vc.items(): mass[b(sz)] += sz
    ax.bar(np.arange(6)+(k-1)*w, mass, w, color=cc, linewidth=0,
           label="%d%% (%d clusters)" % ({"clu70":70,"clu80":80,"clu90":90}[col], len(vc)))
ax.set_xticks(range(6)); ax.set_xticklabels(BL, fontsize=6.5)
ax.set_xlabel("Cluster size"); ax.set_ylabel("Complexes"); ax.set_ylim(0,72)
ax.legend(frameon=False, fontsize=5.8, loc="upper center", ncol=1,
          handlelength=0.9, borderpad=0.2)
ax.set_title("Cluster-size mass", fontsize=8, pad=3)
ax.text(-0.26,1.02,"C", transform=ax.transAxes, **PL)

# --- D  per-complex max identity to training --------------------------
ax = fig.add_subplot(gs[1,0:3])
for p in ORDER:
    v = np.sort(pc.loc[pc.protocol==p,"max_sim"].values)
    ax.step(v, np.arange(1,len(v)+1)/len(v), where="post",
            color=PC[p], lw=1.7 if p in ("tcr70","random") else 1.0,
            label="%s (n=%d)" % (p, len(v)))
ax.axvline(0.90, ls="--", lw=0.8, color="0.4")
ax.text(0.898, 0.50, "0.90", fontsize=6, color="0.4", ha="right",
        rotation=90, va="center")
ax.set_xlabel("Max identity to any training complex")
ax.set_ylabel("Cumulative fraction"); ax.set_xlim(0.52,1.02); ax.set_ylim(0,1.03)
ax.legend(frameon=False, fontsize=6, loc="upper left", handlelength=1.1)
ax.text(0.97,0.30,"random: median %.3f, %.0f%% \u2265 0.90\ntcr70: max %.3f, 0%% \u2265 0.90"
        % (np.median(pc.loc[pc.protocol=="random","max_sim"]),
           100*(pc.loc[pc.protocol=="random","max_sim"]>=0.9).mean(),
           pc.loc[pc.protocol=="tcr70","max_sim"].max()),
        transform=ax.transAxes, fontsize=6, va="top", ha="right", linespacing=1.45)
ax.set_title("Test-set proximity to training", fontsize=8, pad=3)
ax.text(-0.10,1.02,"D", transform=ax.transAxes, **PL)

# --- E  duplicate-family leakage --------------------------------------
ax = fig.add_subplot(gs[1,3:6])
g = au.groupby("protocol").agg(n=("n_test","sum"),
        lc=("leak_conservative","mean")).reindex(ORDER)
x = np.arange(len(ORDER))
ax.bar(x, g.lc, 0.55, color=[PC[p] for p in ORDER], linewidth=0)
for xi,v in zip(x, g.lc):
    ax.text(xi, v+0.16, f"{v:.1f}", ha="center", fontsize=6.5, color="0.25")
ax.set_xticks(x)
ax.set_xticklabels(["%s\n%d/174" % (p, g.n[p]) for p in ORDER],
                   fontsize=6.5, linespacing=1.4)
ax.set_ylabel("Split duplicate families per fold")
ax.set_ylim(0, max(g.lc)*1.25)
ax.set_title("Duplicate-family leakage", fontsize=8, pad=3)
ax.text(0.02,0.90,"complexes evaluated shown below each protocol",
        transform=ax.transAxes, fontsize=5.8, color="0.45")
ax.text(-0.10,1.02,"E", transform=ax.transAxes, **PL)

fig.savefig("figures/fig2_dataset_v3.png", dpi=600, facecolor="white")
fig.savefig("figures/fig2_dataset_v3.pdf", facecolor="white")
print("wrote figures/fig2_dataset_v3.png / .pdf")
print("exact twins (NN >= 0.999): %d / %d = %.1f%%"
      % ((nn>=0.999).sum(), len(nn), 100*(nn>=0.999).mean()))

print("\n--- caption numbers ---")
print("complexes %d | clusters 70/80/90%%: %d/%d/%d"
      % (len(m), m.clu70.nunique(), m.clu80.nunique(), m.clu90.nunique()))
print("all pairs %d | median %.3f | >=0.90 %.2f%%"
      % (len(off), np.median(off), 100*(off>=0.9).mean()))
print("nearest-neighbour | median %.3f | >=0.90 %.1f%% | >=0.95 %.1f%% | min %.3f"
      % (np.median(nn), 100*(nn>=0.9).mean(), 100*(nn>=0.95).mean(), nn.min()))
vc_ = m.dup_family.value_counts(); ind_ = (m.dup_family.map(vc_)>1).values
print("median NN: in dup family %.3f (n=%d) | singleton %.3f (n=%d)"
      % (np.median(nn[ind_]), ind_.sum(), np.median(nn[~ind_]), (~ind_).sum()))
print("TRAV12 %d (%.1f%%) | TRBV6 %d (%.1f%%)"
      % (m.trav.value_counts()['TRAV12'], 100*m.trav.value_counts()['TRAV12']/len(m),
         m.trbv.value_counts()['TRBV6'], 100*m.trbv.value_counts()['TRBV6']/len(m)))
print()
print(pc.groupby("protocol").max_sim.agg(["count","mean","median","max"]).round(3).to_string())
print()
print(g.round(2).to_string())
