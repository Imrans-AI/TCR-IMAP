#!/usr/bin/env python3
"""Step 1 of the p84v2 figure build: statistics before any figure. NO training, inference or re-scoring.

1. Derive the pairwise per-complex arrays for p84v2 by executing scripts/24c_derive_pairwise_summaries_p82v1.py's
   own source on pairmaps_tcr70_p84v2.npz (per-complex AUPRC of the stored out-of-fold predictions, the same
   per_complex function the training run used) -> pairwise_ours_p84v2.npz, percomplex_tcr70_paired_p84v2.tsv.
2. Provenance: every recorded p84v2 value that has per-complex arrays is reproduced at 2,000 resamples / seed 0
   (benchmark, protocol_main, ablations, peptide). Any mismatch -> STOP (exit 1).
3. cluster_macro_p84v2.csv: the equal-cluster-weight estimator (figstats.cluster_macro_bootstrap) is first shown to
   reproduce the recorded cluster_macro_p82v1.csv from the p82v1 arrays, then applied to p84v2 with the same n_boot/seed.
4. 20,000 resamples, seed 42 (figstats.pooled_cluster_bootstrap) for every displayed comparison; Holm across the eight
   ablations and Bonferroni-simultaneous 99.375 % intervals from the same resamples.
5. p82v1 vs p84v2 headline table, both computed by this same code from their own arrays.
Outputs: results/tables/{pairwise_ours,percomplex_tcr70_paired,cluster_macro}_p84v2.*, results/review/figures_p84v2/*.tsv
usage: python scripts/55_stats_p84v2.py
"""
import contextlib
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "."); sys.path.insert(0, "src")
from src.eval.metrics import per_complex  # noqa: E402
from viz.figstats import pooled_cluster_bootstrap as pcb, cluster_macro_bootstrap as cmb, holm, parse_cell  # noqa: E402
from src.data import p84v2 as P  # noqa: E402

T = "results/tables/"; OUT = Path("results/review/figures_p84v2"); OUT.mkdir(parents=True, exist_ok=True)
md = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
PROTOS = ["random", "allele", "tcr70", "double", "temporal"]
COMP = ["Struct-TCS", "Seq-TCS", "PosFreq", "NN"]
bad = []
def same(x, y, tol=5.1e-5): return abs(float(x) - float(y)) <= tol

# ---------------- 1. pairwise per-complex arrays (24c's own code) ----------------
if not Path(T + "pairwise_ours_p84v2.npz").exists():
    src = Path("scripts/24c_derive_pairwise_summaries_p82v1.py").read_text()
    for old, new, n in (("pairmaps_tcr70_p82v1.npz", "pairmaps_tcr70_p84v2.npz", 2), ("pairwise_ours_p82v1.npz", "pairwise_ours_p84v2.npz", 2),
                        ("percomplex_tcr70_paired_p82v1.tsv", "percomplex_tcr70_paired_p84v2.tsv", 2)):
        assert src.count(old) == n, old; src = src.replace(old, new)
    with contextlib.redirect_stdout(io.StringIO()) as buf:
        exec(compile(src, "24c [p84v2]", "exec"), {"__name__": "__main__"})
    print("derived pairwise per-complex arrays (24c on p84v2):", buf.getvalue().strip().splitlines()[-1])


def pairwise(tag, protocol):
    pm = np.load(T + f"pairmaps_{protocol}_{tag}.npz", allow_pickle=True)
    i = np.where(pm["covered"])[0]
    Y, S = pm["Y"].astype(np.float32)[i], pm["S"].astype(np.float32)[i]
    return (per_complex(Y, pm["pred_model"][i], S)["auprc"], per_complex(Y, pm["pred_freq"][i], S)["auprc"],
            pm["clu70"][i], pm["pdb"][i])


# ---------------- 2. provenance at 2,000 / seed 0 ----------------
prov = []
bench = pd.read_csv(T + "benchmark_p84v2.csv").set_index("protocol")
for p in PROTOS:
    z = np.load(T + f"percomplex_{p}_p84v2.npz", allow_pickle=True)
    for k in COMP:
        rec = parse_cell(bench.loc[p, f"vs_{k}"]); r = pcb(z["Ours_auprc"], z[f"{k}_auprc"], z["clu70"], 2000, 0)
        ok = all(same(a, b) for a, b in zip((r["diff"], r["lo"], r["hi"], r["p"]), rec))
        prov.append(dict(table="benchmark_p84v2", item=f"per-residue {p}: Ours vs {k}", recorded=f"{rec}", reproduced=f"{(round(r['diff'],4), round(r['lo'],4), round(r['hi'],4), round(r['p'],4))}", ok=ok))
pmain = pd.read_csv(T + "protocol_main_p84v2.csv").set_index("protocol")
PW = {}
for p in PROTOS:
    a, b, cl, pdbs = pairwise("p84v2", p); PW[p] = (a, b, cl, pdbs)
    r = pcb(a, b, cl, 2000, 0); row = pmain.loc[p]
    ok = all(same(x, y) for x, y in ((np.nanmean(a), row.model_auprc), (np.nanmean(b), row.freq_auprc), (r["diff"], row["diff"]),
                                     (r["lo"], row.ci_lo), (r["hi"], row.ci_hi), (r["p"], row.pval)))
    prov.append(dict(table="protocol_main_p84v2", item=f"pairwise {p}: model vs pair prior", recorded=f"{row['diff']} [{row.ci_lo},{row.ci_hi}] p={row.pval}",
                     reproduced=f"{r['diff']:.4f} [{r['lo']:.4f},{r['hi']:.4f}] p={r['p']:.4f}", ok=ok))
pw = np.load(T + "pairwise_ours_p84v2.npz", allow_pickle=True)
prov.append(dict(table="pairwise_ours_p84v2 (24c)", item="equals the tcr70 per-complex model AUPRC used above", recorded="", reproduced="",
                 ok=bool(np.allclose(pw["auprc"], PW["tcr70"][0], equal_nan=True))))
CACHE = T + "_partial/p84v2/29b_ablations/"
abl = pd.read_csv(T + "ablations_p84v2.csv").set_index("ablation")
full = np.load(CACHE + "full.npz"); cl70 = md.loc[[str(x) for x in pw["pdb"]], "clu70"].values
prov.append(dict(table="29b cache", item="full model per-complex AUPRC equals pairwise_ours_p84v2", recorded="", reproduced="",
                 ok=bool((full["idx"] == np.arange(174)).all() and np.allclose(full["auprc"], pw["auprc"], equal_nan=True))))
ABL = ["no_mhc", "no_tcr", "shuffle_tcr", "shuffle_mhc", "no_germline", "no_cdr3", "no_cross", "dot_head"]
AB = {}
for k in ABL:
    a = np.load(CACHE + f"{k}.npz"); AB[k] = a["auprc"]
    r = pcb(full["auprc"], a["auprc"], cl70, 2000, 0); rr = abl.loc[k]
    ok = all(abs(x - y) < 1e-12 for x, y in ((r["diff"], rr["diff"]), (r["lo"], rr.lo), (r["hi"], rr.hi), (r["p"], rr.p)))
    prov.append(dict(table="ablations_p84v2", item=k, recorded=f"{rr['diff']:+.4f} [{rr.lo:+.4f},{rr.hi:+.4f}] p={rr.p}",
                     reproduced=f"{r['diff']:+.4f} [{r['lo']:+.4f},{r['hi']:+.4f}] p={r['p']}", ok=ok))
# peptide: per-complex from stored predictions, rows as 34b defines them (p84v2 segment layout)
pp = np.load(T + "peptide_preds_p84v2.npz", allow_pickle=True); lab = P.load()
pi = pp["idx"].astype(int); Yp, Sp = lab["Y"][pi], lab["S"][pi]; clp = md.loc[[lab["pdb"][i] for i in pi], "clu70"].values
lens = lab["seg_lens"]; off = 0; gm, cm = [], []
for s in P.SEG_ORDER:
    (gm if s in P.GERMLINE_SEGS else cm).extend(range(off, off + lens[s])); off += lens[s]
PEP = {}
for scope, rows_ in (("all", None), ("germline(CDR1/2)", gm), ("CDR3", cm)):
    for tag in ("no_peptide", "peptide", "peptide_shuf"):
        Ysub, Ssub, Psub = (Yp, Sp, pp[tag][pi]) if rows_ is None else (Yp[:, rows_], Sp[:, rows_], pp[tag][pi][:, rows_])
        PEP[(scope, tag)] = per_complex(Ysub, Psub, Ssub)["auprc"]
pa = pd.read_csv(T + "peptide_aware_p84v2.csv").set_index("scope")
for scope, a_tag, rowname in (("all", "peptide", "peptide"), ("all", "peptide_shuf", "peptide_shuf"),
                              ("germline(CDR1/2)", "peptide", "germline(CDR1/2)"), ("CDR3", "peptide", "CDR3")):
    r = pcb(PEP[(scope, a_tag)], PEP[(scope, "no_peptide")], clp, 2000, 0); row = pa.loc[rowname]
    ok = all(same(x, y) for x, y in ((r["diff"], row["diff"]), (r["lo"], row.lo), (r["hi"], row.hi), (r["p"], row.pval)))
    prov.append(dict(table="peptide_aware_p84v2", item=rowname, recorded=f"{row['diff']} [{row.lo},{row.hi}] p={row.pval}",
                     reproduced=f"{r['diff']:.4f} [{r['lo']:.4f},{r['hi']:.4f}] p={r['p']:.4f}", ok=ok))
te = np.load(T + "teim_res_pairwise_p84v2.npz", allow_pickle=True)
assert list(te["pdb"]) == list(pw["pdb"])
tm = np.load(T + "pairmaps_tcr70_p84v2.npz", allow_pickle=True)
te_re = per_complex(tm["Y"].astype(np.float32), te["pred"], tm["S"].astype(np.float32))["auprc"]
prov.append(dict(table="teim_res_pairwise_p84v2", item="stored per-complex AUPRC equals AUPRC of its stored predictions", recorded=f"{np.nanmean(te['auprc']):.4f}",
                 reproduced=f"{np.nanmean(te_re):.4f}", ok=bool(np.allclose(te["auprc"], te_re, equal_nan=True))))
PROV = pd.DataFrame(prov); PROV.to_csv(OUT / "provenance_2k_seed0.tsv", sep="\t", index=False)
print(f"provenance: {int(PROV.ok.sum())}/{len(PROV)} reproduce")
if not PROV.ok.all():
    print(PROV[~PROV.ok].to_string(index=False)); print("STOP: a recorded p84v2 value does not reproduce"); sys.exit(1)

# ---------------- 3. cluster-macro: verify on p82v1, then p84v2 ----------------
def cm_rows(tag, seed):
    z = np.load(T + f"percomplex_tcr70_{tag}.npz", allow_pickle=True); cl = z["clu70"]
    a, b, clp_, _ = pairwise(tag, "tcr70")
    rows = []
    def macro(x, c): uc, inv = np.unique(c, return_inverse=True); return float((np.bincount(inv, weights=x) / np.bincount(inv)).mean())
    for name, x in (("Ours", z["Ours_auprc"]), ("Struct-TCS", z["Struct-TCS_auprc"]), ("PosFreq", z["PosFreq_auprc"])):
        rows.append(dict(task="per-residue", protocol="tcr70", model=name, n_complexes=len(x), n_clusters=len(np.unique(cl)),
                         complex_macro=round(float(np.nanmean(x)), 4), cluster_macro=round(macro(x, cl), 4)))
    for k in ("Struct-TCS", "PosFreq"):
        r = cmb(z["Ours_auprc"], z[f"{k}_auprc"], cl, 2000, seed)
        rows.append(dict(task="per-residue", protocol="tcr70", model=f"Ours_minus_{k}", n_complexes=r["n"], n_clusters=r["k"],
                         paired_diff_vs_GRIP_complex_macro=round(float(np.nanmean(z["Ours_auprc"] - z[f"{k}_auprc"])), 4),
                         paired_diff_vs_GRIP_cluster_macro=round(r["diff"], 4), ci_lo_cluster_macro=round(r["lo"], 4),
                         ci_hi_cluster_macro=round(r["hi"], 4), p_cluster_macro=round(r["p"], 4), n_boot=2000))
    for name, x in (("model(GRIP)", a), ("freq(prior)", b)):
        rows.append(dict(task="pairwise", protocol="tcr70", model=name, n_complexes=len(x), n_clusters=len(np.unique(clp_)),
                         complex_macro=round(float(np.nanmean(x)), 4), cluster_macro=round(macro(x, clp_), 4)))
    r = cmb(a, b, clp_, 2000, seed)
    rows.append(dict(task="pairwise", protocol="tcr70", model="model_minus_freq", n_complexes=r["n"], n_clusters=r["k"],
                     paired_diff_vs_GRIP_complex_macro=round(float(np.nanmean(a - b)), 4), paired_diff_vs_GRIP_cluster_macro=round(r["diff"], 4),
                     ci_lo_cluster_macro=round(r["lo"], 4), ci_hi_cluster_macro=round(r["hi"], 4), p_cluster_macro=round(r["p"], 4), n_boot=2000))
    return pd.DataFrame(rows)

ref = pd.read_csv(T + "cluster_macro_p82v1.csv")
cols = ["complex_macro", "cluster_macro", "paired_diff_vs_GRIP_complex_macro", "paired_diff_vs_GRIP_cluster_macro",
        "ci_lo_cluster_macro", "ci_hi_cluster_macro", "p_cluster_macro"]
seed_ok = None
for seed in (0, 42, 1):
    mine = cm_rows("p82v1", seed)
    m = ref.merge(mine, on=["task", "protocol", "model"], suffixes=("_ref", ""))
    diffs = [(r.model, c, getattr(r, c + "_ref"), getattr(r, c)) for r in m.itertuples() for c in cols
             if not (pd.isna(getattr(r, c + "_ref")) and pd.isna(getattr(r, c))) and not same(getattr(r, c + "_ref"), getattr(r, c), 5.1e-4)]
    if not diffs and len(m) == len(ref):
        seed_ok = seed; break
print(f"cluster-macro estimator reproduces cluster_macro_p82v1.csv from the p82v1 arrays at 2,000 resamples, seed {seed_ok}: "
      f"{seed_ok is not None}" + ("" if seed_ok is not None else f" (differences: {diffs[:6]})"))
if seed_ok is None:
    print("STOP: cluster-macro estimator does not reproduce the recorded p82v1 values"); sys.exit(1)
CM84 = cm_rows("p84v2", seed_ok)
CM84 = CM84.reindex(columns=list(ref.columns))
if not Path(T + "cluster_macro_p84v2.csv").exists():
    CM84.to_csv(T + "cluster_macro_p84v2.csv", index=False)

# ---------------- 4. 20,000 resamples, seed 42 ----------------
def s20(a, b, cl):
    r = pcb(a, b, cl); return r
st = []
def add(group, item, r, extra=None):
    st.append(dict(group=group, comparison=item, n=r.get("n"), clusters=r.get("k"), diff=r["diff"], lo=r["lo"], hi=r["hi"], p=r["p"], **(extra or {})))
for p in PROTOS:
    z = np.load(T + f"percomplex_{p}_p84v2.npz", allow_pickle=True)
    for k in COMP:
        add("per-residue", f"{p}: TCR-IMAP vs {k}", s20(z["Ours_auprc"], z[f"{k}_auprc"], z["clu70"]))
for p in PROTOS:
    a, b, cl, _ = PW[p]; add("pairwise", f"{p}: TCR-IMAP vs pair-frequency prior", s20(a, b, cl),
                             dict(wins=f"{int((a > b).sum())}/{len(a)}"))
add("pairwise", "tcr70: TCR-IMAP vs TEIM-Res", s20(pw["auprc"], te["auprc"], cl70), dict(wins=f"{int((pw['auprc'] > te['auprc']).sum())}/174"))
zt = np.load(T + "percomplex_tcr70_p84v2.npz", allow_pickle=True)
for k in ("Struct-TCS", "PosFreq"):
    r = cmb(zt["Ours_auprc"], zt[f"{k}_auprc"], zt["clu70"]); add("cluster-macro", f"per-residue tcr70: TCR-IMAP vs {k}", r)
r = cmb(PW["tcr70"][0], PW["tcr70"][1], PW["tcr70"][2]); add("cluster-macro", "pairwise tcr70: TCR-IMAP vs pair-frequency prior", r)
SIM_Q = 0.05 / 8 / 2
abl20 = {}
for k in ABL:
    r = pcb(full["auprc"], AB[k], cl70, return_draws=True)
    lo, hi = np.quantile(r.pop("draws"), [SIM_Q, 1 - SIM_Q]); r["sim_lo"], r["sim_hi"] = float(lo), float(hi); abl20[k] = r
ph = holm([abl20[k]["p"] for k in ABL])
for k, h in zip(ABL, ph):
    r = abl20[k]; add("ablation (full - ablated)", k, r, dict(p_holm=float(h), sim_lo=r["sim_lo"], sim_hi=r["sim_hi"],
                                                           survives_holm=bool(h < 0.05), sim_excludes_zero=not (r["sim_lo"] <= 0 <= r["sim_hi"])))
for scope, a_tag, name in (("all", "peptide", "peptide vs no peptide"), ("all", "peptide_shuf", "shuffled peptide vs no peptide"),
                           ("germline(CDR1/2)", "peptide", "germline rows: peptide vs no peptide"), ("CDR3", "peptide", "CDR3 rows: peptide vs no peptide")):
    add("peptide", name, s20(PEP[(scope, a_tag)], PEP[(scope, "no_peptide")], clp))
ST = pd.DataFrame(st); ST.to_csv(OUT / "stats_20k_seed42.tsv", sep="\t", index=False)

# ---------------- 5. p82v1 vs p84v2 headline table (same code, own arrays) ----------------
H = []
def both(name, f):
    r1, r2 = f("p82v1"), f("p84v2")
    H.append(dict(headline=name, p82v1=r1, p84v2=r2))
def fmt(r): return f"{r['diff']:+.4f} [{r['lo']:+.4f}, {r['hi']:+.4f}] p={r['p']:.4f}"
def pr_cmp(tag, prot, k):
    z = np.load(T + f"percomplex_{prot}_{tag}.npz", allow_pickle=True); return fmt(pcb(z["Ours_auprc"], z[f"{k}_auprc"], z["clu70"]))
def mean_of(tag, prot, k):
    z = np.load(T + f"percomplex_{prot}_{tag}.npz", allow_pickle=True); return f"{np.nanmean(z[f'{k}_auprc']):.4f}"
for prot in PROTOS:
    both(f"per-residue {prot}: TCR-IMAP macro-AUPRC", lambda t, p=prot: mean_of(t, p, "Ours"))
    for k in ("Struct-TCS", "PosFreq"):
        both(f"per-residue {prot}: TCR-IMAP vs {k} (20k)", lambda t, p=prot, k=k: pr_cmp(t, p, k))
for prot in PROTOS:
    both(f"pairwise {prot}: TCR-IMAP vs pair prior (20k)", lambda t, p=prot: fmt(pcb(*pairwise(t, p)[:3])))
both("pairwise tcr70: TCR-IMAP / prior macro-AUPRC", lambda t: "%.4f / %.4f" % tuple(np.nanmean(x) for x in pairwise(t, "tcr70")[:2]))
both("pairwise tcr70: wins over prior", lambda t: (lambda a, b: f"{int((a > b).sum())}/{len(a)}")(*pairwise(t, "tcr70")[:2]))
def teim_cmp(t):
    pwt = np.load(T + f"pairwise_ours_{t}.npz", allow_pickle=True); tt = np.load(T + f"teim_res_pairwise_{t}.npz", allow_pickle=True)
    return fmt(pcb(pwt["auprc"], tt["auprc"], pwt["clu70"])) + f"; TEIM-Res {np.nanmean(tt['auprc']):.4f}; wins {int((pwt['auprc'] > tt['auprc']).sum())}/174"
both("pairwise tcr70: TCR-IMAP vs TEIM-Res (20k)", teim_cmp)
def cmacro(t):
    a, b, cl, _ = pairwise(t, "tcr70"); return fmt(cmb(a, b, cl))
both("cluster-macro pairwise tcr70 vs prior (20k)", cmacro)
def ablate(t):
    c = T + f"_partial/{t}/29b_ablations/"; f0 = np.load(c + "full.npz")["auprc"]
    pwt = np.load(T + f"pairwise_ours_{t}.npz", allow_pickle=True); clt = md.loc[[str(x) for x in pwt["pdb"]], "clu70"].values
    ps = [pcb(f0, np.load(c + f"{k}.npz")["auprc"], clt)["p"] for k in ABL]; hh = holm(ps)
    return ", ".join(k for k, h in zip(ABL, hh) if h < 0.05) or "none"
both("ablations surviving Holm across 8 (20k)", ablate)
SIMC = pd.read_csv(T + "similarity_curve_p82v1.csv").set_index("threshold"), pd.read_csv(T + "similarity_curve_p84v2.csv").set_index("threshold")
for thr in (1.01, 0.9, 0.7):
    H.append(dict(headline=f"identity sweep < {thr if thr < 1 else 'none'} (recorded 2k)",
                  p82v1=f"{SIMC[0].loc[thr, 'diff']:+.4f} [{SIMC[0].loc[thr, 'lo']:+.4f}, {SIMC[0].loc[thr, 'hi']:+.4f}] p={SIMC[0].loc[thr, 'p']:.3f}",
                  p84v2=f"{SIMC[1].loc[thr, 'diff']:+.4f} [{SIMC[1].loc[thr, 'lo']:+.4f}, {SIMC[1].loc[thr, 'hi']:+.4f}] p={SIMC[1].loc[thr, 'p']:.3f}"))
HD = pd.DataFrame(H); HD.to_csv(OUT / "headline_p82v1_vs_p84v2.tsv", sep="\t", index=False)
with pd.option_context("display.width", 250, "display.max_colwidth", 80, "display.float_format", "{:.4f}".format):
    print("\n20k / seed 42:"); print(ST.to_string(index=False))
    print("\nheadline p82v1 vs p84v2:"); print(HD.to_string(index=False))
    print("\ncluster_macro_p84v2.csv:"); print(CM84.to_string(index=False))
