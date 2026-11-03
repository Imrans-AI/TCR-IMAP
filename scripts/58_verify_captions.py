#!/usr/bin/env python3
"""Verify every number in docs/manuscript/FIGURE_CAPTIONS.md against the committed tables.

For each figure the script computes, from the tables the figure itself reads, every value its caption quotes, formatted
as the caption prints it (sign, U+2212 minus, decimals, thousands separator). It then extracts every number from the
caption and fails unless each one is either such a table-derived value or a declared method constant for that figure
(a cutoff, an IMGT range, a resampling count, a threshold); the constants are listed per figure below, with what they
are. Numbers inside identifiers (CDR1, TRAV12, 5yxn, 1061A, α1, S3) are not numbers and are skipped.

Reads committed tables (and, where a figure already recomputes from them, the stored per-complex arrays) only;
nothing is trained or scored. usage: python scripts/58_verify_captions.py
"""
import re, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from viz.figstats import pooled_cluster_bootstrap, cluster_macro_bootstrap, parse_cell

T, R = "results/tables/", "results/review/"
RD = dict(float_precision="round_trip")
MI = "−"
def s4(x): return f"{x:+.4f}".replace("-", MI)
def s3(x): return f"{x:+.3f}".replace("-", MI)
def s2(x): return f"{x:+.2f}".replace("-", MI)
def u(x, nd): return f"{x:.{nd}f}".replace("-", MI)
def p3(p): return "0.001" if p <= 0.0005 else f"{p:.3f}"      # pstr: floored values print as "p < 0.001"
def ci4(d): return [s4(d["diff"]), s4(d["lo"]), s4(d["hi"])]
def n(x): return f"{int(x):,}"

V = {}          # figure -> {token: source}
C = {}          # figure -> {token: what the constant is}
def add(fig, tokens, src):
    for t in (tokens if isinstance(tokens, (list, tuple)) else [tokens]):
        V.setdefault(fig, {})[t] = src
def const(fig, tokens, what):
    for t in tokens.split():
        C.setdefault(fig, {})[t] = what
ALL = ["Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5", "Figure 6", "Figure 7"] + [f"Figure S{i}" for i in range(1, 10)]

md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t")
N = len(md); K = md.clu70.nunique()
assert (N, K) == (174, 69)
S20 = pd.read_csv(R + "figures_p84v2/stats_20k_seed42.tsv", sep="\t", **RD)
def s20(grp, comp):
    r = S20[(S20.group == grp) & (S20.comparison == comp)]
    assert len(r) == 1, (grp, comp); return r.iloc[0]
for f in ALL:
    add(f, [n(N)], "complex_metadata.tsv: 174 complexes")
    add(f, [str(K)], "complex_metadata.tsv: clu70 clusters")
    const(f, "95", "95 % confidence interval")

# ======================================================================= Figure 2
# ======================================================================= Figure 1
# Figure 1 is the author-approved two-panel figure (figures/main/Figure_1.SOURCE.md), installed on
# 2 November 2026. Its caption quotes only the contact cutoff, the two input dimensions and the
# sigmoid output range; the dataset counts, protocol table and predeclaration chronology that the
# superseded five-panel caption carried now live in Methods, whose numbers Figure 2 to Figure 7 and
# the SI tables already pin to the same artifacts. The row and column counts are still read from the
# frozen tensor so a change in the representation would fail this check.
F = "Figure 1"
_z = np.load("data/processed/pair_maps_p84v2.npz", allow_pickle=True)
assert _z["Y"].shape == (N, 75, 84)
add(F, [str(_z["Y"].shape[1]), str(_z["Y"].shape[2])],
    "pair_maps_p84v2.npz: TCR rows and MHC columns")
const(F, "4.5", "contact cutoff (A)")
const(F, "0 1", "sigmoid output range (contact probability)")

F = "Figure 2"
vc_a, vc_b = md.trav.value_counts(), md.trbv.value_counts()
add(F, [str(vc_a.iloc[0]), str(vc_b.iloc[0])], "complex_metadata.tsv trav/trbv value_counts (TRAV12, TRBV6)")
assert vc_a.index[0] == "TRAV12" and vc_b.index[0] == "TRBV6"
Sim = np.load(T + "tcr_identity_matrix.npy")
iu = np.triu_indices(N, 1)
add(F, n(len(iu[0])), "tcr_identity_matrix.npy: complex pairs")
add(F, u(100 * (Sim[iu] >= 0.90).mean(), 1), "tcr_identity_matrix.npy: pairs >= 0.90 (%)")
nn = np.where(np.eye(N, dtype=bool), -1, Sim).max(1)
add(F, u(100 * (nn >= 0.90).mean(), 1), "tcr_identity_matrix.npy: nearest neighbor >= 0.90 (%)")
tw = int((nn >= 0.999).sum()); add(F, [str(tw), u(100 * tw / N, 1)], "tcr_identity_matrix.npy: exact twins")
add(F, [str(md[c].nunique()) for c in ("clu70", "clu80", "clu90")], "complex_metadata.tsv clu70/80/90")
mx = pd.read_csv(T + "percomplex_maxsim.tsv", sep="\t", **RD)
rnd = mx[mx.protocol == "random"].max_sim
add(F, [u(rnd.median(), 3), u(100 * (rnd >= 0.90).mean(), 1)], "percomplex_maxsim.tsv random: median, % >= 0.90")
t70 = mx[mx.protocol == "tcr70"].max_sim
add(F, u(t70.max(), 3), "percomplex_maxsim.tsv tcr70 max"); assert (t70 < 0.90).all()
lk = pd.read_csv(T + "leakage_audit_dual.tsv", sep="\t", **RD)
ev = {p: len(np.load(T + f"percomplex_{p}_p84v2.npz", allow_pickle=True)["pdb"]) for p in ("random", "double", "allele", "temporal", "tcr70")}
for p, nd in (("random", 1), ("double", 2), ("allele", 2), ("temporal", 0), ("tcr70", 0)):
    sub = lk[lk.protocol == p]
    sub = sub[sub.n_test >= 3]                                  # the evaluation's minimum fold size (allele fold of 2)
    add(F, u(sub.leak_conservative.mean(), nd) if nd else str(int(round(sub.leak_conservative.mean()))),
        f"leakage_audit_dual.tsv {p}: mean leak_conservative over evaluated folds")
    add(F, str(ev[p]), f"percomplex_{p}_p84v2.npz: complexes evaluated")
add(F, str(int(md.allele.isna().sum())), "complex_metadata.tsv: allele missing")
const(F, "70 80 90", "clustering thresholds (%)")
const(F, "0.90", "similarity threshold drawn in B and D")
const(F, "27 38 56 65 105 124 118", "IMGT position ranges of the slot sequences")

# ======================================================================= Figure 3
F = "Figure 3"
sc = pd.read_csv(T + "similarity_curve_p84v2.csv", **RD).set_index("threshold")
a, b = sc.loc[1.01], sc.loc[0.7]
add(F, [u(a.model - b.model, 4), u(a.model, 4), u(b.model, 4), u(a.freq - b.freq, 4), u(a.freq, 4), u(b.freq, 4)],
    "similarity_curve_p84v2.csv: model/freq at none and 0.70")
assert u(a.model - b.model, 4) == "0.3042" and u(a.freq - b.freq, 4) == "0.0566"
add(F, [s4(a["diff"]), p3(a.p), s4(b["diff"]), s4(b.lo), s4(b.hi), p3(b.p)],
    "similarity_curve_p84v2.csv (superseded 2,000-resample values, retained in the caption as the contrast)")
# authoritative sweep statistics: 20,000 resamples, seed 42 (results/review/sweep_authoritative/)
_auth = pd.read_csv(T + "similarity_curve_p84v2_authoritative.csv", **RD).set_index("threshold")
for _t in _auth.index:
    _r = _auth.loc[_t]
    add(F, [s4(_r["diff"]), s4(_r.lo), s4(_r.hi), p3(_r.p)],
        f"similarity_curve_p84v2_authoritative.csv threshold {_t}: diff/lo/hi/p at 20,000 seed 42")
assert all(np.diff(sc.sort_index(ascending=False)["diff"].values) < 0), "advantage decreases at every step"
ang = pd.DataFrame([json.loads(l) for l in open(T + "angle_runs_p84v2.jsonl")])
for s in ("random", "allele", "trbv"):
    d = ang[ang.split_by == s]
    add(F, [u(d.mae.min(), 2), u(d.mae.max(), 2), u(d.pearson.min(), 2), u(d.pearson.max(), 2),
            s2(d.r2.min()), s2(d.r2.max())], f"angle_runs_p84v2.jsonl {s}: MAE, Pearson r, R2 ranges")
add(F, [u(ang.baseline_mae.iloc[0], 1), str(int(ang.n.iloc[0]))], "angle_runs_p84v2.jsonl: median baseline, n")
const(F, "0.70", "strictest threshold"); const(F, "2,000", "recorded bootstrap resamples")
const(F, "0", "R² = 0 reference (mean predictor)")
# panel E: the training-size-matched control (results/review/control_estimands/)
ce = pd.read_csv(R + "control_estimands/control_estimands.tsv", sep="\t", **RD).set_index("threshold")
sd = pd.read_csv(R + "control_estimands/similarity_diagnostics.tsv", sep="\t", **RD)
for thr in (0.95, 0.85, 0.70):
    r = ce.loc[thr]
    add(F, [s4(r.principal_mean_of_draw_differences), s4(r.principal_lo), s4(r.principal_hi),
            s4(r.filtered_delta), s4(r.contrast_filtered_minus_sizematched),
            s4(r.contrast_lo), s4(r.contrast_hi)],
        f"control_estimands.tsv thr {thr}: principal, filtered, paired contrast")
    f = sd[(sd.regime == "similarity-filtered") & (np.isclose(sd.threshold, thr))].iloc[0]
    ts = [int(x) for x in str(f.train_sizes_per_fold).split(",")]
    add(F, [str(min(ts)), str(max(ts))], f"similarity_diagnostics.tsv thr {thr}: retained training sizes per fold")
unr = sd[sd.regime == "unrestricted"].iloc[0]
uts = [int(x) for x in str(unr.train_sizes_per_fold).split(",")]
add(F, [str(min(uts)), str(max(uts))], "similarity_diagnostics.tsv: unrestricted training sizes per fold")
add(F, s4(sc.loc[1.01]["diff"]), "similarity_curve_p84v2.csv: unrestricted advantage (panel E reference line)")
const(F, "0.95 0.85", "the two further thresholds with a size-matched control")
const(F, "20,000 42", "control resamples; seed"); const(F, "0.001", "pstr floor")
const(F, "5", "fixed size-matched draws")

# ======================================================================= Figure 4
F = "Figure 4"
bench = pd.read_csv(T + "benchmark_p84v2.csv", **RD).set_index("protocol")
for p in bench.index:
    add(F, [str(int(bench.loc[p, "n"])), u(bench.loc[p, "Ours"], 3), u(bench.loc[p, "Struct-TCS"], 3)],
        f"benchmark_p84v2.csv {p}: n, Ours, Struct-TCS")
    r = s20("per-residue", f"{p}: TCR-IMAP vs Struct-TCS"); add(F, [s4(r["diff"]), p3(r.p)], f"stats_20k {p} vs Struct-TCS")
r = s20("per-residue", "tcr70: TCR-IMAP vs Struct-TCS"); add(F, ci4(r), "stats_20k tcr70 vs Struct-TCS CI")
for comp, w in (("tcr70: TCR-IMAP vs pair-frequency prior", "112/174"), ("tcr70: TCR-IMAP vs TEIM-Res", "143/174")):
    r = s20("pairwise", comp); add(F, ci4(r) + [p3(r.p)], f"stats_20k pairwise {comp}")
    h, t = map(int, r.wins.split("/")); assert f"{h}/{t}" == w
    add(F, [str(h), u(100 * h / t, 1)], f"stats_20k wins {comp}")
for comp in ("per-residue tcr70: TCR-IMAP vs Struct-TCS", "per-residue tcr70: TCR-IMAP vs PosFreq",
             "pairwise tcr70: TCR-IMAP vs pair-frequency prior"):
    add(F, p3(s20("cluster-macro", comp).p), f"stats_20k cluster-macro {comp}")
add(F, p3(s20("per-residue", "tcr70: TCR-IMAP vs PosFreq").p), "stats_20k tcr70 vs PosFreq")
F20 = pd.read_csv(R + "validation_fix/F_structcs_differences_20k.tsv", sep="\t", **RD)
for _, r in F20[(F20.comparison == "TCR-IMAP vs Struct-TCS") & (F20.config == "cluster-disjoint (reported)")].iterrows():
    add(F, ci4(r) + [p3(r.p), str(int(r.n_test))], f"F_structcs_differences_20k {r.cutoff} {r.setting}")
st = pd.read_csv(T + "structcs_setting_p84v2_6A.csv").set_index("setting")
add(F, str(int(st.n_shared_with_comparator.iloc[0])), "structcs_setting_p84v2_6A.csv n_shared_with_comparator")
const(F, "4.5 6", "the comparator's contact cutoffs (Å)"); const(F, "20,000 42", "resamples, seed")

# ======================================================================= Figure 5
F = "Figure 5"
z = np.load(T + "contact_freq_gdomain182_descriptive_p84v2.npz", allow_pickle=True)
lab = [str(x) for x in z["tcr_position_labels"]]
add(F, str(len(lab)), "contact_freq npz: TCR rows")
ins = [l for l in lab if re.search(r"\d[A-Z]$", l)]
add(F, [str(len(lab) - len(ins)), str(len(ins))], "contact_freq npz: occupied IMGT positions / CDR3 insertion rows")
reg = pd.read_csv(T + "gdomain182_region_classification_p84v2.csv")
add(F, str(len(z["mhc_positions"])), "contact_freq npz: G-domain columns")
win = [str(x) for x in z["mhc_positions"] if 50 <= int(re.match(r"\d+", str(x))[0]) <= 90 or 1050 <= int(re.match(r"\d+", str(x))[0]) <= 1090]
a1 = [w for w in win if int(re.match(r"\d+", w)[0]) < 1000]
add(F, [str(len(win)), str(len(a1)), str(len(win) - len(a1))], "contact_freq npz: window columns (all, α1, α2)")
Y, S = z["Y"], z["S"]
res_rows = (S.sum(axis=(2)) > 0).sum(0)                         # complexes in which each TCR row is resolved
rare = res_rows < 5
assert sorted(lab[i] for i in np.where(rare)[0]) == ["cdr3b:111B", "cdr3b:112B", "cdr3b:112C"], [lab[i] for i in np.where(rare)[0]]
colw = np.array([str(x) in win for x in z["mhc_positions"]])
Fw = z["F"][:, colw]; Fw = np.where(rare[:, None], -1, Fw)
i, j = np.unravel_index(np.nanargmax(Fw), Fw.shape)
assert lab[i] == "cdr1b:36" and np.array(win)[j] == "76"
add(F, ["36", "76", u(Fw[i, j], 3)], "contact_freq npz: hottest pair over rows resolved in >= 5 complexes")
hot = (np.where(rare[:, None], 0, z["F"][:, colw])).sum(0); hot = hot / hot.max()
hp = [w for w, h in zip(win, hot) if h > 0.30]
assert hp == ["62", "65", "66", "68", "69", "72", "73", "76", "1061A", "1062", "1065", "1066", "1069", "1073"], hp
add(F, [str(len(hp))] + [h for h in hp if h.isdigit()], "contact_freq npz: positions > 0.30 of maximum")
cg = int(Y.sum()); cw = int(Y[:, :, colw].sum())
add(F, [n(cw), n(cg), u(100 * cw / cg, 2)], "contact_freq npz: window vs G-domain contacts")
rg = dict(zip(reg.imgt_position.astype(str), reg.region))
fl = np.array([rg[str(x)] == "beta-sheet floor" for x in z["mhc_positions"]])
assert len(fl) == Y.shape[2]
add(F, n(int(S[:, :, fl].sum())), "contact_freq npz + region classification: resolved pairs on the β-sheet floor")
assert int(Y[:, :, fl].sum()) == 0
gd = pd.read_csv(T + "germline_diagonal_p84v2.csv", **RD).set_index("loop")
add(F, n(int(gd.loc["cdr2b", "resolved_a2"])), "germline_diagonal_p84v2.csv cdr2b resolved_a2"); assert gd.loc["cdr2b", "contacts_a2"] == 0
add(F, [p3(gd.loc["cdr3b", "p_bootstrap"]), u(gd.p_holm.max(), 3)], "germline_diagonal_p84v2.csv p_bootstrap (CDR3β), p_holm")
assert (gd.p_holm == 0.003).all() and (gd.drop("cdr3b").p_bootstrap <= 0.0005).all()
assert dict(gd.prefers) == {"cdr1a": "alpha2", "cdr2a": "alpha2", "cdr3a": "alpha1", "cdr1b": "alpha1", "cdr2b": "alpha1", "cdr3b": "alpha2"}
const(F, "4.5", "contact cutoff (Å)"); const(F, "50 90 1050 1090 1000", "window definition (IMGT)")
const(F, "0.30", "stick threshold"); const(F, "0 2,000 42", "β-floor count stated as zero; resamples; seed")
const(F, "0.001", "pstr floor")

# ======================================================================= Figure 6
F = "Figure 6"
ab = pd.read_csv(T + "ablations_p84v2.csv", **RD)
add(F, u(float(ab[ab.ablation == "full"].auprc.iloc[0]), 4), "ablations_p84v2.csv full auprc")
add(F, u(100 * (1 - 0.05 / 8), 3), "Bonferroni simultaneous level for 8 ablations")
for k in ("shuffle_tcr", "no_germline", "shuffle_mhc", "no_tcr", "no_mhc", "no_cross", "dot_head"):
    r = s20("ablation (full - ablated)", k); add(F, [s4(r["diff"]), p3(r.p_holm)], f"stats_20k ablation {k}")
r = s20("pairwise", "tcr70: TCR-IMAP vs pair-frequency prior")
add(F, ci4(r) + [p3(r.p), "112", u(100 * 112 / 174, 1)], "stats_20k pairwise vs prior")
pr = pd.read_csv(T + "percomplex_tcr70_paired_p84v2.tsv", sep="\t", **RD).set_index("pdb")
assert pr.delta.idxmax() == "3qfj" and pr.delta.idxmin() == "3qeq"
sites = pd.read_csv("figures/pymol/fig6_sites_p84v2.tsv", sep="\t")
for p in ("3qfj", "3qeq"):
    add(F, [u(pr.loc[p, "grip"], 3), u(pr.loc[p, "prior"], 3), s4(pr.loc[p, "delta"]), str(int(pr.loc[p, "npos"]))],
        f"percomplex_tcr70_paired_p84v2.tsv {p}")
    o, q = sites[(sites.pdb == p) & (sites.set == "observed")], sites[(sites.pdb == p) & (sites.set == "predicted")]
    assert len(o) == len(q) == int(pr.loc[p, "npos"])
    add(F, [str(int(q.observed.sum())), str(len(set(o.mhc_imgt) & set(q.mhc_imgt))), str(len(set(o.mhc_imgt)))],
        f"fig6_sites_p84v2.tsv {p}: predicted pairs observed, MHC sites shared / observed")
const(F, "0.05", "Holm significance level"); const(F, "4.5 84", "contact cutoff (Å); window columns")
const(F, "20,000 42 0.001", "resamples; seed; pstr floor")

# ======================================================================= Figure 7
F = "Figure 7"
rho = pd.read_csv(T + "md_occupancy_rho.csv", **RD); strata = pd.read_csv(T + "md_occupancy_strata.csv", **RD)
ts = pd.read_csv(T + "md_characterisation_timeseries.csv", **RD); summ = pd.read_csv(T + "md_characterisation_summary.csv", **RD)
man = pd.read_csv(T + "md_occupancy_manifest.csv").set_index("item")["value"]
assert man["permutation"].startswith("within-system, 10,000") and "two-sided" in man["permutation"]
nfr = summ.n_window_frames.unique(); assert list(nfr) == [9001], nfr
assert list(ts.dropna(subset=["interface_rmsd"]).groupby(["pdb", "rep"]).size().unique()) == [9001]
add(F, ["9,001", "10,000"], "md_characterisation_summary.csv n_window_frames; md_occupancy_manifest permutation")
pool = rho[rho.level == "pooled"].set_index("score")
add(F, [n(pool.loc["TCR-IMAP", "n_crystal_negative"]), str(int(pool.loc["TCR-IMAP", "n_nonzero_occupancy"])),
        n(pool.loc["TCR-IMAP", "n_crystal_negative"] - pool.loc["TCR-IMAP", "n_nonzero_occupancy"]),
        u(pool.loc["pair-frequency prior", "rho"], 2), u(pool.loc["TCR-IMAP", "rho"], 2)], "md_occupancy_rho.csv pooled")
sysr = rho[rho.level.isin(["system", "replicate"])].pivot_table(index=["pdb", "rep"], columns="score", values="rho")
assert (sysr["pair-frequency prior"] > sysr["TCR-IMAP"]).all(), "prior > TCR-IMAP in every system and replicate"
sp = strata[strata.level == "pooled"].set_index(["stratum", "score"])
for s_ in ("4.5-6 A", "6-8 A", "> 8 A"):
    add(F, [n(sp.loc[(s_, "TCR-IMAP"), "n"]), u(sp.loc[(s_, "TCR-IMAP"), "rho"], 2), u(sp.loc[(s_, "pair-frequency prior"), "rho"], 2),
            str(int(sp.loc[(s_, "TCR-IMAP"), "n_nonzero_occupancy"]))], f"md_occupancy_strata.csv pooled {s_}")
assert int(sp.loc[("> 8 A", "TCR-IMAP"), "n_occupancy_ge_025"]) == 0
assert rho.interpretable.all() and strata[strata.level != "replicate"].interpretable.all()
assert min(rho.n_nonzero_occupancy.min(), strata[strata.level != "replicate"].n_nonzero_occupancy.min()) >= 10
ss = pd.read_csv(T + "md_score_score_spearman.csv", **RD)
add(F, [u(ss[ss.level == "system"].spearman_rho.min(), 2), u(ss[ss.level == "system"].spearman_rho.max(), 2),
        u(ss[ss.level == "pooled"].spearman_rho.iloc[0], 2)], "md_score_score_spearman.csv")
tk = pd.read_csv(T + "md_occupancy_topk.csv", **RD)
enr = tk.frac_top_k / tk.base_rate
add(F, [str(tk.k.min()), str(tk.k.max()), str(int(np.floor(enr.min()))) if False else str(int(round(enr.min()))),
        str(int(round(enr.max()))), str(tk.n_top_k_ge.min()), str(tk.n_top_k_ge.max())], "md_occupancy_topk.csv k, enrichment, counts")
dk = pd.read_csv(T + "md_occupancy_docked_frames.csv", **RD)
d5 = dk[dk.pdb == "5yxn"].set_index("rep").docked_fraction_count * 100
assert (d5 < 50).sum() == 2
add(F, [f"{d5.loc[1]:.0f}", f"{d5.loc[3]:.0f}"], "md_occupancy_docked_frames.csv 5yxn rep 1, 3 (count)")
c1 = pd.read_csv(T + "md_occupancy_check1_positive.csv", **RD)
sm = c1[c1.rep.astype(str) == "system"].set_index("pdb")["mean"]
assert sm.idxmin() == "5yxn"; add(F, u(sm.loc["5yxn"], 2), "md_occupancy_check1_positive.csv 5yxn system mean")
for p in ("5yxn", "7n1f", "3vxs", "8eo8"):
    add(F, u(summ[summ.pdb == p].backbone_rmsd_mean.mean(), 2), f"md_characterisation_summary.csv {p} backbone RMSD")
bb = ts.backbone_rmsd
add(F, [u(np.floor(bb.min() * 10) / 10, 1), u(np.ceil(bb.max() * 10) / 10, 1)], "md_characterisation_timeseries.csv backbone RMSD range (rounded outward)")
const(F, "200 20 4.5 6 8 0.25 50 0.05", "run length, window start (ns); distance cutoffs and strata (Å); occupancy threshold; docked 50 %; fill threshold")
const(F, "10", "interpretable minimum of non-zero occupancies")

# ======================================================================= Figure S1
F = "Figure S1"
pa = pd.read_csv(T + "peptide_aware_p84v2.csv", **RD)
for _, r in pa.iterrows():
    add(F, [s4(r["diff"]), s4(r.lo), s4(r.hi), p3(r.pval)], f"peptide_aware_p84v2.csv {r.scope}")
assert ((pa.lo < 0) & (pa.hi > 0)).all()
const(F, "2,000", "recorded resamples")

# ======================================================================= Figure S2
F = "Figure S2"
for p in bench.index:
    add(F, str(int(bench.loc[p, "n"])), f"benchmark_p84v2.csv {p} n")
const(F, "1.5", "whisker length (× IQR)")

# ======================================================================= Figure S3
F = "Figure S3"
add(F, "9,001", "md_characterisation_timeseries.csv frames per replicate")
sasa_n = ts.dropna(subset=["buried_sasa"]).groupby(["pdb", "rep"]).size().unique(); assert list(sasa_n) == [181]
add(F, "181", "md_characterisation_timeseries.csv SASA points per replicate")
s5 = summ[summ.pdb == "5yxn"].set_index("rep")
assert summ.interface_rmsd_mean.idxmax() == s5.index.map(lambda r: summ[(summ.pdb == "5yxn") & (summ.rep == r)].index[0])[0]
add(F, u(s5.loc[1, "interface_rmsd_mean"], 2), "md_characterisation_summary.csv 5yxn rep1 interface RMSD mean (highest of 12)")
ir1 = ts[(ts.pdb == "5yxn") & (ts.rep == 1)].interface_rmsd.dropna()
assert (ir1 > 8).mean() > 0.25, (ir1 > 8).mean()
add(F, ["8", u(ir1.max(), 1)], "md_characterisation_timeseries.csv 5yxn rep1 interface RMSD: > 8 A in > 25 % of frames; maximum")
t5 = ts[(ts.pdb == "5yxn") & (ts.rep == 2)].dropna(subset=["buried_sasa"]).assign(t=lambda d: d.time_ps / 1000)
first, last = t5[t5.t <= 30].buried_sasa.mean(), t5[t5.t >= 190].buried_sasa.mean()     # 20-30 and 190-200 ns, inclusive (CHARACTERISATION_FIX_REPORT)
print(f"5yxn rep2 SASA first/last 10 ns: {first:.1f} / {last:.1f}")
add(F, [n(round(first)), n(round(last))], "timeseries 5yxn rep2 SASA, means of 20-30 and 190-200 ns")
sys_sasa = summ.groupby("pdb").buried_sasa_mean.mean(); assert sys_sasa.idxmin() == "5yxn"
add(F, n(round(sys_sasa.loc["5yxn"])), "md_characterisation_summary.csv 5yxn mean buried SASA")
add(F, [u(s5.loc[r, "rmsf_cdr3b"], 2) for r in (1, 2, 3)], "md_characterisation_summary.csv 5yxn rmsf_cdr3b")
const(F, "1 2 3", "replicate numbers")
const(F, "200 20 27 38 56 65 105 117 1.4 1 3.5 150 10", "run length, window, frame spacing (ns/ps); IMGT loop ranges; probe (Å); SASA stride (ns); H-bond criteria; 10 ns window")

# ======================================================================= Figure S4
F = "Figure S4"
mac = pd.read_csv(R + "validation_fix/B3_macro_auprc.tsv", sep="\t", **RD)
add(F, [s4(c) for c in mac[mac.method == "TCR-IMAP"].change], "B3_macro_auprc.tsv TCR-IMAP change")
kd = pd.read_csv(R + "validation_fix/B3_key_differences_20k.tsv", sep="\t", **RD)
for _, r in kd[kd.protocol == "tcr70"].iterrows():
    add(F, ci4(r) + [p3(r.p)], f"B3_key_differences_20k tcr70 {r.comparison} {r.split}")
lg = open(R + "validation_fix/queue_console_randomsplit_deterministic.log").read()
cnt = [int(x) for x in re.findall(r"\((\d+) clu70 clusters on both sides\)", lg)]
assert len(cnt) == 18
add(F, ["18", str(min(cnt)), str(max(cnt)), f"{np.median(cnt):g}"], "queue_console_randomsplit_deterministic.log straddle counts")
for p in bench.index:
    add(F, str(int(bench.loc[p, "n"])), f"benchmark_p84v2.csv {p} n")
const(F, "15 0", "random validation fraction (%); straddling clusters under the cluster-disjoint split (asserted in B2_split_dryrun.log)")
assert "cluster-disjoint in every one: True" in open(R + "validation_fix/B2_split_dryrun.log").read()
const(F, "20,000 42 0.001", "resamples; seed; pstr floor")

# ======================================================================= Figure S5
F = "Figure S5"
for g, c in (("per-residue", "tcr70: TCR-IMAP vs Struct-TCS"), ("per-residue", "tcr70: TCR-IMAP vs PosFreq"),
             ("pairwise", "tcr70: TCR-IMAP vs pair-frequency prior"), ("cluster-macro", "per-residue tcr70: TCR-IMAP vs Struct-TCS"),
             ("cluster-macro", "per-residue tcr70: TCR-IMAP vs PosFreq"), ("cluster-macro", "pairwise tcr70: TCR-IMAP vs pair-frequency prior")):
    r = s20(g, c); add(F, ci4(r) + [p3(r.p)], f"stats_20k {g} {c}")
const(F, "20,000 42 0.001", "resamples; seed; pstr floor")

# ======================================================================= Figure S6
F = "Figure S6"
for _, r in F20[F20.comparison == "TCR-IMAP vs Struct-TCS"].iterrows():
    add(F, ci4(r) + [p3(r.p), str(int(r.n_test))], f"F_structcs_differences_20k {r.cutoff} {r.setting} {r.config}")
stab = pd.read_csv(R + "pre_md_run/B_perseed_stability.tsv", sep="\t", **RD).set_index(["cutoff", "setting"])
assert not stab.sign_flip.any()
add(F, [s3(stab.loc[("6A", "their_random"), f"seed{s}"]) for s in (42, 43, 44)], "B_perseed_stability.tsv 6A their_random seeds")
add(F, "88", "structcs_setting_p84v2_6A.csv n_shared_with_comparator")
const(F, "4.5 6", "the comparator's contact cutoffs (Å)"); const(F, "20,000 42 0.001", "resamples; seed; pstr floor")

# ======================================================================= Figure S7
F = "Figure S7"
pm = pd.read_csv(T + "protocol_main.csv", **RD).set_index("protocol").loc["tcr70"]
bv = pd.read_csv(T + "benchmark_v2.csv", **RD).set_index("protocol").loc["tcr70"]
sc28 = pd.read_csv(T + "similarity_curve.csv", **RD).set_index("threshold").loc[0.7]
pa28 = pd.read_csv(T + "peptide_aware.csv", **RD).set_index("scope").loc["CDR3"]
add(F, [s4(pm["diff"]), p3(pm.pval)], "protocol_main.csv tcr70")
for c in ("vs_Struct-TCS", "vs_PosFreq"):
    d_, lo_, hi_, p_ = parse_cell(bv[c]); add(F, [s4(d_), p3(p_)], f"benchmark_v2.csv tcr70 {c}")
add(F, [s4(sc28["diff"]), p3(sc28.p)], "similarity_curve.csv 0.70"); add(F, [s4(pa28["diff"]), p3(pa28.pval)], "peptide_aware.csv CDR3")
for g, c in (("pairwise", "tcr70: TCR-IMAP vs pair-frequency prior"), ("per-residue", "tcr70: TCR-IMAP vs Struct-TCS"),
             ("per-residue", "tcr70: TCR-IMAP vs PosFreq")):
    r = s20(g, c); add(F, [s4(r["diff"]), p3(r.p)], f"stats_20k {g} {c}")
add(F, [s4(sc.loc[0.7, "diff"]), p3(sc.loc[0.7, "p"])], "similarity_curve_p84v2.csv 0.70")
pa84 = pa.set_index("scope").loc["CDR3"]; add(F, [s4(pa84["diff"]), p3(pa84.pval)], "peptide_aware_p84v2.csv CDR3")
L = np.load("data/processed/pair_maps.npz", allow_pickle=True)["Y"]; add(F, str(L.shape[2]), "pair_maps.npz MHC positions")
add(F, "84", "pair_maps_p84v2.npz MHC columns"); assert np.load("data/processed/pair_maps_p84v2.npz", allow_pickle=True)["Y"].shape[2] == 84
const(F, "0.70 2,000 20,000 42 0.001", "sweep threshold; resamples; seed; pstr floor")

# ======================================================================= Figure S8
F = "Figure S8"
pr8 = pd.read_csv(T + "percomplex_tcr70_paired_p84v2.tsv", sep="\t")
add(F, str(int(pr8.allele.isna().sum())), "percomplex_tcr70_paired_p84v2.tsv allele missing")
va, vb = pr8.allele.dropna().value_counts(), pr8.trbv.fillna("NA").value_counts()
add(F, [str(int((va < 3).sum())), str(int(va[va < 3].sum())), str(int((vb < 3).sum())), str(int(vb[vb < 3].sum()))],
    "percomplex_tcr70_paired_p84v2.tsv groups with n < 3 (alleles; TRBV incl. missing)")
assert vb.get("NA", 0) == 1

# ======================================================================= Figure S9
F = "Figure S9"
co = pd.read_csv(T + "similarity_curve_composition_p84v2.csv", **RD).set_index("threshold")
fo = pd.read_csv(T + "similarity_curve_composition_folds_p84v2.csv", **RD)
add(F, [str(fo.n_test.min()), str(fo.n_test.max())], "composition folds: test per fold")
add(F, [str(int(co.loc[1.01, "n_train_min"])), str(int(co.loc[1.01, "n_train_max"])), str(int(co.loc[0.7, "n_train_min"])),
        str(int(co.loc[0.7, "n_train_max"])), u(co.loc[0.7, "achieved_max_sim_max"], 3), u(co.loc[0.7, "achieved_max_sim_mean"], 3)],
    "similarity_curve_composition_p84v2.csv")
add(F, [s4(sc.loc[0.7, "diff"]), p3(sc.loc[0.7, "p"])], "similarity_curve_p84v2.csv 0.70")
const(F, "0.70 2,000", "strictest threshold; recorded resamples")

# ======================================================================= check the captions
text = open("docs/manuscript/FIGURE_CAPTIONS.md").read()
secs = dict(re.findall(r"^## (Figure S?\d+)\n(.*?)(?=^## |\Z)", text, flags=re.S | re.M))
assert set(ALL) <= set(secs), set(ALL) - set(secs)
TOK = re.compile(r"(?<![\w./])[−+]?\d[\d,]*(?:\.\d+)?(?![\w/²])")
bad_total, report = 0, []
for f in ALL:
    body = re.sub(r"Figures? S?\d+[A-F]?(?:/[A-F])?(?:,? ?[A-F](?![\w]))*", " ", secs[f])     # figure references
    body = body.replace("R²", "R2X")                                                      # 'R² = 0' keeps its 0
    toks = [t.rstrip(",") for t in TOK.findall(body)]
    ok = {**V.get(f, {}), **{k: "constant: " + v for k, v in C.get(f, {}).items()}}
    bad = sorted({t for t in toks if t not in ok})
    used = sorted({t for t in toks if t in V.get(f, {})})
    report.append(f"{f}: {len(toks)} numbers, {len(set(toks))} distinct; {len(used)} traced to tables, "
                  f"{len({t for t in toks if t in C.get(f, {}) and t not in V.get(f, {})})} declared constants; "
                  f"UNVERIFIED: {bad if bad else 'none'}")
    for t in sorted(set(toks)):
        src = V.get(f, {}).get(t) or ("constant: " + C[f][t] if t in C.get(f, {}) else "UNVERIFIED")
        report.append(f"    {t:>10}  {src}")
    bad_total += len(bad)
print("\n".join(report))
print(f"\n{'ALL CAPTION NUMBERS VERIFIED' if bad_total == 0 else f'{bad_total} UNVERIFIED NUMBER(S)'}")
sys.exit(1 if bad_total else 0)
