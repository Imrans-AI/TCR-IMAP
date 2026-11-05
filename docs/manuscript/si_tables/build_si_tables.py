#!/usr/bin/env python3
"""Build the Supporting Information tables S1-S8.

Reads committed tables and reports, and the stored per-complex arrays (verified by SHA-256 where a checksum is
committed). No training, inference, scoring or GPU. The only statistics computed here are paired cluster-bootstrap
differences recomputed from stored per-complex arrays with the existing figstats.pooled_cluster_bootstrap (20,000
resamples, seed 42); every one is asserted equal to the committed value that already records it, and the build stops
on any mismatch. Parameter counts are taken by instantiating each model class on the CPU (no forward pass).

Writes only new files under docs/manuscript/si_tables/: S1-S8 as Markdown + CSV, and INDEX.md.
usage (repo root): CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=4 nice -n 15 python docs/manuscript/si_tables/build_si_tables.py
"""
import ast
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "src"))
from viz.figstats import pooled_cluster_bootstrap as pcb, holm  # noqa: E402

OUT = ROOT / "docs/manuscript/si_tables"
T = ROOT / "results/tables"
RV = ROOT / "results/review"
for f in OUT.glob("S*_*"):
    raise SystemExit(f"refusing to overwrite existing {f}")

MINUS = "−"

# Display name for the partition whose frozen key is "double". It is an intersection-selected
# partition, not a double holdout: its held-out groups remain in training
# (results/review/V3_PART_A_VERIFICATION.md A1). The frozen key is NEVER renamed -- protocols.json
# and every stored artifact keep "double" -- only what the tables print changes.
PLAB = {"double": "intersection (legacy)"}
def plab(p):
    return PLAB.get(p, p)

PENDING = []


def sg(x, nd=4):
    return f"{x:+.{nd}f}".replace("-", MINUS)


def fx(x, nd=4):
    return f"{x:.{nd}f}".replace("-", MINUS)


def ci(lo, hi, nd=4):
    return f"[{sg(lo, nd)}, {sg(hi, nd)}]"


def pstr(p):
    """Project convention (src/viz/style.pstr): '< 0.001' at or below the 2,000-resample floor, else 3 dp."""
    return "< 0.001" if p <= 0.0005 else f"{p:.3f}"


def n_(x):
    return f"{int(x):,}"


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def same(a, b, tol=1e-12):
    return abs(float(a) - float(b)) <= tol


def check(cond, msg):
    if not cond:
        raise SystemExit(f"CHECK FAILED: {msg}")
    print(f"  ok  {msg}")


def write(stem, md_text, csvs):
    (OUT / f"{stem}.md").write_text(md_text.rstrip() + "\n")
    for name, df in csvs.items():
        df.to_csv(OUT / f"{name}.csv", index=False)
    print(f"wrote {stem}.md + {', '.join(n + '.csv' for n in csvs)}")


# ---------------------------------------------------------------- inputs and their integrity
print("== inputs")
md = pd.read_csv(ROOT / "data/processed/complex_metadata.tsv", sep="\t")
check(sha(ROOT / "data/processed/complex_metadata.tsv").startswith("c14641f75b25dc264b2655791e66d2e5"),
      "complex_metadata.tsv sha256 = frozen value in CORRECTION_PLAN_P84V2.md")
check(sha(ROOT / "data/splits_v2/protocols.json").startswith("3e48a094e94d9758c54187ea924e23ca"),
      "protocols.json sha256 = frozen value in CORRECTION_PLAN_P84V2.md")
sums = {l.split()[1]: l.split()[0] for l in (RV / "validation_fix/clusterdisjoint_run1/SHA256SUMS").read_text().splitlines() if l.strip()}
PROTOS = ["random", "double", "allele", "temporal", "tcr70"]
for p in PROTOS:
    check(sha(T / f"percomplex_{p}_p84v2.npz") == sums[f"percomplex_{p}_p84v2.npz"],
          f"percomplex_{p}_p84v2.npz = committed SHA256SUMS (validation_fix/clusterdisjoint_run1, the reported run)")
ST = pd.read_csv(RV / "figures_p84v2/stats_20k_seed42.tsv", sep="\t")


def strow(group, comparison):
    r = ST[(ST.group == group) & (ST.comparison == comparison)]
    assert len(r) == 1, (group, comparison)
    return r.iloc[0]


# ================================================================= S1 benchmark
print("== S1")
bench = pd.read_csv(T / "benchmark_p84v2.csv").set_index("protocol")
METH = [("Ours", "TCR-IMAP"), ("Struct-TCS", "Struct-TCS"), ("Seq-TCS", "Seq-TCS"),
        ("PosFreq", "Position-frequency prior"), ("NN", "Nearest-neighbor transfer")]
STNAME = {"Struct-TCS": "Struct-TCS", "Seq-TCS": "Seq-TCS", "PosFreq": "PosFreq", "NN": "NN"}
s1a, s1a_csv = [], []
for p in PROTOS:
    z = np.load(T / f"percomplex_{p}_p84v2.npz", allow_pickle=True)
    n, k = len(z["pdb"]), len(np.unique(z["clu70"]))
    for key, label in METH:
        apr, auc = z[f"{key}_auprc"], z[f"{key}_auc"]
        m_apr = float(np.nanmean(apr)); m_auc = float(np.nanmean(auc))
        check(round(m_apr, 4) == round(float(bench.loc[p, key]), 4), f"S1 per-residue {p} {label} macro-AUPRC = benchmark_p84v2.csv")
        check(int(np.isfinite(apr).sum()) == n and int(np.isfinite(auc).sum()) == n, f"S1 per-residue {p} {label}: AUPRC and AUC defined for all {n} complexes")
        row = dict(task="per-residue", protocol=p, method=label, n_complexes=n, n_clusters=k,
                   n_auprc_defined=int(np.isfinite(apr).sum()), macro_auprc=m_apr, macro_auc=m_auc,
                   diff_vs_tcrimap=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan, resamples="")
        cell = ""
        if key != "Ours":
            r = pcb(z["Ours_auprc"], apr, z["clu70"])
            s = strow("per-residue", f"{p}: TCR-IMAP vs {STNAME[key]}")
            check(same(r["diff"], s["diff"]) and same(r["lo"], s["lo"]) and same(r["hi"], s["hi"]) and r["p"] == s["p"],
                  f"S1 per-residue {p} TCR-IMAP vs {label} at 20k = stats_20k_seed42.tsv")
            row.update(diff_vs_tcrimap=r["diff"], ci_lo=r["lo"], ci_hi=r["hi"], p=r["p"], resamples="20,000 (seed 42)")
            cell = f"{sg(r['diff'])} {ci(r['lo'], r['hi'])}", pstr(r["p"])
        s1a_csv.append(row)
        s1a.append([plab(p) if key == "Ours" else "", n if key == "Ours" else "", k if key == "Ours" else "", label,
                    fx(m_apr), fx(m_auc), *(cell if cell else ("n/a", "n/a"))])

pm = pd.read_csv(T / "protocol_main_p84v2.csv").set_index("protocol")
pw = np.load(T / "pairwise_ours_p84v2.npz", allow_pickle=True)
te = np.load(T / "teim_res_pairwise_p84v2.npz", allow_pickle=True)
pt = pd.read_csv(T / "percomplex_tcr70_paired_p84v2.tsv", sep="\t")
check(list(pt.pdb) == list(pw["pdb"]) == list(te["pdb"]), "S1 pairwise tcr70 arrays share complex order")
check(np.allclose(pt.grip.values, pw["auprc"], atol=1e-12, equal_nan=True), "S1 percomplex_tcr70_paired grip = pairwise_ours auprc")
# Pairwise ROC AUC for every protocol, scored from the stored prediction arrays in
# results/tables/pairmaps_{protocol}_p84v2.npz (pred_model, pred_freq, Y, S, covered). These were
# previously left as [PENDING] on the belief that no per-complex pairwise AUC was stored; the
# per-complex AUC is not stored, but the PREDICTIONS are, so the AUC is computed here exactly as
# the macro-AUPRC in the same table is. Validated below against the two values already committed.
from eval.metrics import per_complex as _pc   # noqa: E402
_LAB = np.load(ROOT / "data/processed/pair_maps_p84v2.npz", allow_pickle=True)
PWAUC = {}
for _p in PROTOS:
    _pm = np.load(T / f"pairmaps_{_p}_p84v2.npz", allow_pickle=True)
    check([str(x) for x in _pm["pdb"]] == list(md.pdb), f"S1 pairmaps {_p} identifier order = metadata")
    check(bool(np.array_equal(_pm["clu70"], md.clu70.values)), f"S1 pairmaps {_p} clu70 = metadata")
    check(bool(np.array_equal(_pm["Y"].astype(np.int8), _LAB["Y"])), f"S1 pairmaps {_p} labels = frozen build")
    check(bool(np.array_equal(_pm["S"].astype(np.int8), _LAB["S"])), f"S1 pairmaps {_p} masks = frozen build")
    _i = np.where(_pm["covered"])[0]
    _Y, _S = _pm["Y"].astype(np.float32)[_i], _pm["S"].astype(np.float32)[_i]
    _mm, _ff = _pc(_Y, _pm["pred_model"][_i], _S), _pc(_Y, _pm["pred_freq"][_i], _S)
    PWAUC[_p] = (float(np.nanmean(_mm["auc"])), float(np.nanmean(_ff["auc"])))
    check(round(float(np.nanmean(_mm["auprc"])), 4) == round(float(pm.loc[_p, "model_auprc"]), 4),
          f"S1 pairwise {_p} macro-AUPRC rescored from predictions = protocol_main_p84v2.csv")
check(round(PWAUC["tcr70"][0], 4) == round(float(np.nanmean(pw["auc"])), 4),
      "S1 pairwise tcr70 AUC rescored from predictions = stored pairwise_ours_p84v2.npz")

s1b, s1b_csv = [], []
for p in PROTOS:
    n = int(pm.loc[p, "n"]); s = strow("pairwise", f"{p}: TCR-IMAP vs pair-frequency prior")
    check(int(s["n"]) == n, f"S1 pairwise {p} n = protocol_main_p84v2.csv")
    auc_m = fx(PWAUC[p][0]); auc_f = fx(PWAUC[p][1])
    if p == "tcr70":
        r = pcb(pt.grip.values, pt.prior.values, pt.clu70.values)
        check(same(r["diff"], s["diff"]) and same(r["lo"], s["lo"]) and same(r["hi"], s["hi"]) and r["p"] == s["p"],
              "S1 pairwise tcr70 TCR-IMAP vs prior at 20k recomputed from percomplex_tcr70_paired_p84v2.tsv = stats_20k")
        check(round(float(np.nanmean(pw["auprc"])), 4) == round(float(pm.loc[p, "model_auprc"]), 4), "S1 pairwise tcr70 mean = protocol_main")
        check(same(PWAUC[p][0], float(np.nanmean(pw["auc"]))), "S1 pairwise tcr70 AUC = stored array")
    s1b_csv += [dict(task="pairwise", protocol=p, method="TCR-IMAP", n_complexes=n, n_clusters=int(s["clusters"]),
                     macro_auprc=float(pm.loc[p, "model_auprc"]), macro_auc=PWAUC[p][0],
                     diff_vs_tcrimap=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan, wins="", resamples=""),
                dict(task="pairwise", protocol=p, method="Pair-frequency prior", n_complexes=n, n_clusters=int(s["clusters"]),
                     macro_auprc=float(pm.loc[p, "freq_auprc"]), macro_auc=PWAUC[p][1], diff_vs_tcrimap=float(s["diff"]),
                     ci_lo=float(s["lo"]), ci_hi=float(s["hi"]), p=float(s["p"]), wins=s["wins"], resamples="20,000 (seed 42)")]
    s1b.append([plab(p), n, int(s["clusters"]), "TCR-IMAP", fx(pm.loc[p, "model_auprc"]), auc_m, "n/a", "n/a", "n/a"])
    s1b.append(["", "", "", "Pair-frequency prior", fx(pm.loc[p, "freq_auprc"]), auc_f,
                f"{sg(s['diff'])} {ci(s['lo'], s['hi'])}", pstr(s["p"]), s["wins"]])
    if p == "tcr70":
        r = pcb(pw["auprc"], te["auprc"], pw["clu70"]); s2 = strow("pairwise", "tcr70: TCR-IMAP vs TEIM-Res")
        check(same(r["diff"], s2["diff"]) and same(r["lo"], s2["lo"]) and same(r["hi"], s2["hi"]) and r["p"] == s2["p"],
              "S1 pairwise tcr70 TCR-IMAP vs TEIM-Res at 20k = stats_20k")
        wins = f"{int((pw['auprc'] > te['auprc']).sum())}/174"
        check(wins == s2["wins"], "S1 TEIM-Res wins = stats_20k")
        s1b.append(["", "", "", "TEIM-Res", fx(np.nanmean(te["auprc"])), fx(np.nanmean(te["auc"])),
                    f"{sg(r['diff'])} {ci(r['lo'], r['hi'])}", pstr(r["p"]), wins])
        s1b_csv.append(dict(task="pairwise", protocol=p, method="TEIM-Res", n_complexes=174, n_clusters=69,
                            macro_auprc=float(np.nanmean(te["auprc"])), macro_auc=float(np.nanmean(te["auc"])),
                            diff_vs_tcrimap=r["diff"], ci_lo=r["lo"], ci_hi=r["hi"], p=r["p"], wins=wins, resamples="20,000 (seed 42)"))
PENDING[:] = sorted(set(PENDING))

FS = pd.read_csv(RV / "validation_fix/F_structcs_differences_20k.tsv", sep="\t")
s1c, s1c_csv = [], []
for cut, tag, d in (("4.5 Å", "4p5A", "31b_structcs_45A"), ("6 Å", "6A", "31b_structcs_6A")):
    sc = pd.read_csv(T / f"structcs_setting_p84v2_{tag}.csv").set_index("setting")
    for setting, lab in (("their_random", "single random split"), ("our_tcr70", "tcr70")):
        z = np.load(T / f"_partial/p84v2/{d}/{tag}_{setting}.npz", allow_pickle=True)
        pc = json.loads(str(z["__percomplex__"]))
        o = np.array(pc["Ours_auprc"], float); cl = np.array(pc["clu70"])
        for comp in ("Struct-TCS", "Seq-TCS"):
            c = np.array(pc[f"{comp}_auprc"], float)
            check(round(float(c.mean()), 4) == round(float(sc.loc[setting, f"{comp}_auprc"]), 4), f"S1 subset {tag} {setting} {comp} mean = structcs_setting csv")
            r = pcb(o, c, cl)
            f = FS[(FS.cutoff == tag) & (FS.setting == setting) & (FS.comparison == f"TCR-IMAP vs {comp}") & (FS.config == "cluster-disjoint (reported)")].iloc[0]
            check(same(r["diff"], f["diff"]) and same(r["lo"], f["lo"]) and same(r["hi"], f["hi"]) and r["p"] == f["p"],
                  f"S1 subset {tag} {setting} vs {comp} at 20k = F_structcs_differences_20k.tsv")
            s1c_csv.append(dict(cutoff=cut, setting=lab, n_complexes=len(o), n_clusters=len(np.unique(cl)), comparator=comp,
                                tcrimap_macro_auprc=float(o.mean()), comparator_macro_auprc=float(c.mean()),
                                diff=r["diff"], ci_lo=r["lo"], ci_hi=r["hi"], p=r["p"], resamples="20,000 (seed 42)"))
            s1c.append([cut, lab, len(o), len(np.unique(cl)), comp, fx(o.mean()), fx(c.mean()), f"{sg(r['diff'])} {ci(r['lo'], r['hi'])}", pstr(r["p"])])

ang = pd.read_json(T / "angle_runs_p84v2.jsonl", lines=True)
s1d = [[r.split_by, r.seed, r.n, f"{r.mae:.2f}", f"{r.baseline_mae:.2f}", f"{r.pearson:.3f}", f"{r.spearman:.3f}", sg(r.r2, 3)]
       for r in ang.itertuples()]
check(set(ang.n) == {172}, "S1 AngleNet n = 172 in every run")

s1_md = f"""# Table S1. Full benchmark: every method under every protocol

**Table S1.** Per-complex macro-averaged AUPRC and ROC AUC of every method under every partitioning protocol, with the
number of test complexes and of TCR clusters (clu70) they span, and the paired difference from TCR-IMAP. (A) Per-residue
task. (B) Pairwise task. (C) The 88 structures shared with the comparator's dataset. (D) Docking-angle regression
(AngleNet, a separate regressor; not TCR-IMAP).

## A. Per-residue task (cluster-disjoint validation; the reported configuration)

{md_table(["Protocol", "n complexes", "n clusters", "Method", "Macro-AUPRC", "Macro-AUC", "TCR-IMAP − method [95 % CI]", "p"], s1a)}

## B. Pairwise task

{md_table(["Protocol", "n complexes", "n clusters", "Method", "Macro-AUPRC", "Macro-AUC", "TCR-IMAP − method [95 % CI]", "p", "TCR-IMAP higher"], s1b)}

## C. Comparator subset: 88 structures shared with the Struct-TCS dataset (per-residue; cluster-disjoint validation)

{md_table(["Cutoff", "Setting", "n complexes", "n clusters", "Comparator", "TCR-IMAP macro-AUPRC", "Comparator macro-AUPRC", "TCR-IMAP − comparator [95 % CI]", "p"], s1c)}

## D. Docking-angle regression (AngleNet; one row per partition and training seed)

{md_table(["Partition", "Seed", "n", "MAE (°)", "Median-predictor MAE (°)", "Pearson r", "Spearman ρ", "R²"], s1d)}

**Notes.**
- Macro-AUPRC and macro-AUC are means over test complexes of the per-complex value; a complex with no positive or no
  negative resolved candidate has no defined AUPRC or AUC and is omitted from that mean (none is omitted in any cell of A).
- Differences are paired: TCR-IMAP − method, per complex, then averaged. 95 % CI and two-sided p from the pooled clu70
  cluster bootstrap (`src/viz/figstats.pooled_cluster_bootstrap`), **20,000 resamples, seed 42, recomputed by
  `build_si_tables.py` from the stored per-complex arrays** and asserted equal to the committed values; p is floored at
  1/20,000 and reported as < 0.001 at or below 0.0005 (project convention). The p values in the
  `vs_*` columns of `benchmark_p84v2.csv` and `structcs_setting_p84v2_*.csv` are the run's own 2,000-resample values and
  differ slightly; the 20,000-resample values here are the ones reported in the paper.
- "TCR-IMAP higher": complexes on which TCR-IMAP's AUPRC exceeds the comparator's.
- n clusters: distinct clu70 clusters among the evaluated complexes (the bootstrap unit).
- Pairwise ROC AUC is computed here by scoring the stored prediction arrays in
  `results/tables/pairmaps_{{protocol}}_p84v2.npz`, exactly as the macro-AUPRC in the same rows is. Those arrays were
  verified against the frozen labels, masks, identifiers and cluster assignments before use, and the recomputed
  tcr70 values reproduce the committed `pairwise_ours_p84v2.npz` AUC and the macro-AUPRC of every protocol.
- D: R² is relative to the mean predictor; the MAE baseline is the median predictor. AngleNet uses its own five-fold
  random, allele-holdout and TRBV-holdout partitions (not the frozen protocols); 172 of 174 complexes (3vxu and 3w0w have
  no docking angle).

**Sources.** A: `results/tables/percomplex_{{protocol}}_p84v2.npz` (SHA-256 verified against
`results/review/validation_fix/clusterdisjoint_run1/SHA256SUMS`), means checked against `results/tables/benchmark_p84v2.csv`,
differences checked against `results/review/figures_p84v2/stats_20k_seed42.tsv`. B: `results/tables/protocol_main_p84v2.csv`
(means), `stats_20k_seed42.tsv` (differences, wins; recomputed for tcr70 from `percomplex_tcr70_paired_p84v2.tsv`),
`pairwise_ours_p84v2.npz` and `teim_res_pairwise_p84v2.npz` (tcr70 AUC, TEIM-Res). The pairwise differences under
random, intersection (legacy), allele and temporal are taken from `stats_20k_seed42.tsv` as committed (their per-complex arrays were
derived in memory by `scripts/55_stats_p84v2.py` and not persisted). C: `results/tables/structcs_setting_p84v2_{{4p5A,6A}}.csv`,
`results/tables/_partial/p84v2/31b_structcs_{{45A,6A}}/*.npz`, checked against
`results/review/validation_fix/F_structcs_differences_20k.tsv` (cluster-disjoint rows). D: `results/tables/angle_runs_p84v2.jsonl`.
"""
write("S1_benchmark", s1_md, {"S1_benchmark_per_residue_pairwise": pd.DataFrame(s1a_csv + s1b_csv),
                               "S1_benchmark_comparator_subset": pd.DataFrame(s1c_csv), "S1_benchmark_anglenet": ang})

# ================================================================= S2 ablations
print("== S2")
ABL = ["no_mhc", "no_tcr", "shuffle_tcr", "shuffle_mhc", "no_cdr3", "no_germline", "no_cross", "dot_head"]
LAB = {"no_mhc": "MHC residue identity masked", "no_tcr": "TCR residue identity masked",
       "shuffle_tcr": "TCR residues shuffled within each loop", "shuffle_mhc": "MHC residues shuffled across positions",
       "no_cdr3": "CDR3 + FR4 tokens removed", "no_germline": "Germline CDR1/CDR2 tokens removed",
       "no_cross": "Cross-attention removed", "dot_head": "Bilinear head replaced by dot product"}
abl_csv = pd.read_csv(T / "ablations_p84v2.csv").set_index("ablation")
C = T / "_partial/p84v2/29b_ablations"
full = np.load(C / "full.npz", allow_pickle=True)
check(np.allclose(full["auprc"], pw["auprc"], atol=1e-12, equal_nan=True) and list(full["idx"]) == list(range(174)),
      "S2 29b 'full' cache = pairwise_ours_p84v2 (same model, same complexes)")
cl70 = pw["clu70"]
SIM_Q = 0.05 / 8 / 2
res = {}
for k in ABL:
    a = np.load(C / f"{k}.npz", allow_pickle=True)
    check(list(a["idx"]) == list(full["idx"]), f"S2 {k} cache covers the same complexes")
    check(round(float(np.nanmean(a["auprc"])), 6) == round(float(abl_csv.loc[k, "auprc"]), 6), f"S2 {k} macro-AUPRC = ablations_p84v2.csv")
    r = pcb(full["auprc"], a["auprc"], cl70, return_draws=True)
    lo, hi = np.quantile(r.pop("draws"), [SIM_Q, 1 - SIM_Q]); r.update(sim_lo=float(lo), sim_hi=float(hi), auprc=float(np.nanmean(a["auprc"])))
    res[k] = r
ph = holm([res[k]["p"] for k in ABL])
s2, s2_csv = [], []
for k, h in zip(ABL, ph):
    r = res[k]; s = strow("ablation (full - ablated)", k)
    check(all(same(r[x], s[y]) for x, y in (("diff", "diff"), ("lo", "lo"), ("hi", "hi"), ("sim_lo", "sim_lo"), ("sim_hi", "sim_hi")))
          and r["p"] == s["p"] and same(h, s["p_holm"]), f"S2 {k}: diff, CI, simultaneous interval, p and Holm p = stats_20k_seed42.tsv")
    surv = bool(h < 0.05); simx = not (r["sim_lo"] <= 0 <= r["sim_hi"])
    check(surv == simx, f"S2 {k}: simultaneous interval excludes zero iff Holm p < 0.05")
    s2.append([LAB[k], f"`{k}`", fx(r["auprc"]), sg(r["diff"]), ci(r["lo"], r["hi"]), ci(r["sim_lo"], r["sim_hi"]),
               pstr(r["p"]), pstr(h), "yes" if surv else "no"])
    s2_csv.append(dict(ablation=k, description=LAB[k], ablated_macro_auprc=r["auprc"], diff_full_minus_ablated=r["diff"],
                       ci95_lo=r["lo"], ci95_hi=r["hi"], sim99375_lo=r["sim_lo"], sim99375_hi=r["sim_hi"], p_raw=r["p"],
                       p_holm=float(h), survives_holm=surv, n_complexes=r["n"], n_clusters=r["k"], resamples=20000, seed=42))
s2_md = f"""# Table S2. Ablations of TCR-IMAP (tcr70, pairwise task)

**Table S2.** Change in pairwise macro-AUPRC under tcr70 when one input or architectural component is removed, for all
eight ablations. Δ = full model − ablated model (positive, the ablation lowered performance). Full TCR-IMAP:
macro-AUPRC {fx(float(np.nanmean(full['auprc'])))}, n = 174 complexes, 69 clusters.

{md_table(["Ablation", "Code", "Ablated macro-AUPRC", "Δ", "95 % CI", "99.375 % simultaneous interval", "p (raw)", "p (Holm, 8)", "Survives Holm"], s2)}

**Notes.**
- Each ablated model was trained and evaluated exactly as the full model (same folds, epochs, optimizer and seeds 42–44,
  predictions averaged), on the same 174 complexes.
- 95 % CI, simultaneous interval and raw p come from the same 20,000 cluster-bootstrap resamples (pooled clu70, seed 42),
  **recomputed by `build_si_tables.py` from the stored per-complex arrays** and asserted equal to the committed values.
  The simultaneous interval is the Bonferroni 1 − 0.05/8 = 99.375 % percentile interval (quantiles 0.3125 % and 99.6875 %).
- Holm step-down adjustment across the eight ablations. "Survives Holm": Holm-adjusted p < 0.05; in every row this
  coincides with the simultaneous interval excluding zero (asserted).
- p is floored at 1/20,000 and reported as < 0.001 at or below 0.0005. The `p` column of `ablations_p84v2.csv` is the
  run's own 2,000-resample value.
- CDR3 + FR4: the CDR3 token block spans IMGT 105–124, including the J-region framework positions 118–124.

**Sources.** `results/tables/_partial/p84v2/29b_ablations/{{full,{','.join(ABL)}}}.npz` (per-complex AUPRC); clusters
from `results/tables/pairwise_ours_p84v2.npz`; means checked against `results/tables/ablations_p84v2.csv`; all
statistics checked against `results/review/figures_p84v2/stats_20k_seed42.tsv`. Ablation definitions:
`scripts/29b_ablations_p82v1.py` (`apply_ablation`, `NoCross`, `DotHead`), executed on the p84v2 inputs by
`scripts/53_ablations_p84v2.py`.
"""
write("S2_ablations", s2_md, {"S2_ablations": pd.DataFrame(s2_csv)})

# ================================================================= S3 hyperparameters
print("== S3")
import torch  # noqa: E402  (CPU only; classes are instantiated, never called)
torch.set_num_threads(4)
import src.models.grip as g  # noqa: E402
SEG84 = {"cdr1a": 7, "cdr2a": 8, "cdr3a": 22, "cdr1b": 6, "cdr2b": 7, "cdr3b": 25}


def classes(path, names, extra=None):
    src = (ROOT / path).read_text(); tree = ast.parse(src)
    code = "\n".join(ast.get_source_segment(src, n) for n in tree.body if isinstance(n, ast.ClassDef) and n.name in names)
    ns = dict(torch=torch, nn=torch.nn, np=np, VOCAB=g.VOCAB, SEG_ORDER=g.SEG_ORDER, SEG_LENS=SEG84, TCR_LEN=75)
    ns.update(extra or {}); exec(code, ns); return ns


def npar(m):
    return sum(p.numel() for p in m.parameters())


ns19 = classes("scripts/19_train_pairmap.py", {"PairNet"})
ns21 = classes("scripts/21_struct_tcs.py", {"GatedGAT", "AttnBlock", "GATBlock", "OursPerRes", "TCS"})
ns36 = classes("scripts/36b_teim_res_p82v1.py", {"SeqEncoder", "ResBlock2D", "TEIMRes"})
ns16 = classes("scripts/16b_docking_angle_p82v1.py", {"AngleNet"})
PC = {"TCR-IMAP (pairwise)": npar(ns19["PairNet"](mhc_len=84)),
      "TCR-IMAP (per-residue)": npar(ns21["OursPerRes"](mhc_len=84)),
      "Struct-TCS": npar(ns21["TCS"](mode="struct", mhc_len=84)),
      "Seq-TCS": npar(ns21["TCS"](mode="seq", mhc_len=84)),
      "TEIM-Res": npar(ns36["TEIMRes"](84)),
      "AngleNet": npar(ns16["AngleNet"](mhc_len=84))}
for k, v in (("TCR-IMAP (pairwise)", 529310), ("TCR-IMAP (per-residue)", 291073), ("Struct-TCS", 336905), ("Seq-TCS", 303105)):
    check(PC[k] == v, f"S3 parameter count {k} = {v:,} (recorded value)")
print("  parameter counts:", PC)
H = ["Item", "TCR-IMAP (pairwise)", "TCR-IMAP (per-residue)", "Struct-TCS", "Seq-TCS", "TEIM-Res", "AngleNet"]
R = [
    ["Task", "pairwise contact map (75 × 84)", "per-residue (75 rows)", "per-residue", "per-residue", "pairwise contact map", "docking-angle regression"],
    ["Source class (script)", "`PairNet` (19)", "`OursPerRes` (21)", "`TCS(mode='struct')` (21)", "`TCS(mode='seq')` (21)", "`TEIMRes` (36b)", "`AngleNet` (16b)"],
    ["p84v2 runner", "`53_run_all_p84v2.py`", "`53_benchmark_p84v2.py`", "`53_benchmark_p84v2.py`", "`53_benchmark_p84v2.py`", "`53_teim_res_p84v2.py`", "`53_docking_angle_p84v2.py`"],
    ["Inputs", "TCR (75) + MHC (84) tokens", "TCR + MHC tokens", "TCR + MHC tokens + TCR residue graph (Cβ–Cβ ≤ 8 Å; Cα for Gly)", "TCR + MHC tokens", "TCR + MHC tokens", "TCR + MHC tokens; FR4 rows (IMGT 118–124) padded"],
    ["Width d", "128", "96", "64", "64", "64 (embedding); 32 (2D channels)", "96"],
    ["TCR encoder", "2-layer pre-norm Transformer, FF 256", "2-layer pre-norm Transformer, FF 192", "2-layer LSTM → 2 gated multi-head GAT blocks", "2-layer LSTM → 2 self-attention blocks", "2 × Conv1d (k = 3) + BatchNorm + ReLU", "2-layer pre-norm Transformer, FF 192"],
    ["MHC encoder", "1-layer pre-norm Transformer, FF 256", "1-layer pre-norm Transformer, FF 192", "2-layer LSTM → 2 self-attention blocks", "2-layer LSTM → 2 self-attention blocks", "2 × Conv1d (k = 3) + BatchNorm + ReLU", "embedding + positional only"],
    ["Attention heads", "8 (encoders and cross-attention)", "6", "4 (GAT, self- and cross-attention)", "4", "n/a", "6"],
    ["Fusion", "TCR→MHC cross-attention, residual + LayerNorm", "same as pairwise", "cross-attention (TCR query) + FFN", "same as Struct-TCS", "element-wise outer product → 1 × 1 conv → 3 dilated 2D residual blocks (dilation 1, 2, 4)", "TCR→MHC cross-attention, residual + LayerNorm"],
    ["Output head", "bilinear s·⟨W_t h_i, W_m g_j⟩/√d + P_ij + b (P 75 × 84, init 0; s init 1; b init −4)", "per row: Linear–GELU–Dropout–Linear", "per row: Linear", "per row: Linear", "1 × 1 conv", "pooled [TCR mean, MHC mean]: Linear–GELU–Dropout–Linear"],
    ["Dropout", "0.2", "0.25", "0.1 (self-attention blocks)", "0.1", "none", "0.3"],
    ["Loss", "BCE, positive-class weight = neg/pos (resolved training pairs), masked to resolved pairs", "same, per row", "same, per row", "same, per row", "same as pairwise", "smooth L1 on the standardized angle"],
    ["Optimizer", "AdamW", "Adam", "Adam", "Adam", "AdamW", "AdamW"],
    ["Learning rate", *["1 × 10⁻³"] * 6],
    ["Weight decay", "0.01", "1 × 10⁻⁴", "1 × 10⁻⁴", "1 × 10⁻⁴", "0.01", "0.05"],
    ["Gradient clipping (norm)", "1.0", "1.0", "1.0", "1.0", "1.0", "1.0"],
    ["Batch", "full batch", "full batch", "full batch", "full batch", "full batch", "full batch"],
    ["Epochs", "400 (fixed)", "≤ 300", "≤ 300", "≤ 300", "400 (fixed)", "400 (fixed)"],
    ["Validation / early stopping", "none", "cluster-disjoint validation (whole clu70 clusters, ≥ 15 % of training complexes); AUPRC every 5 epochs; stop after 10 evaluations without improvement; best weights restored", "same as per-residue TCR-IMAP", "same", "none", "none"],
    ["Seeds", "42, 43, 44; predictions averaged", "42, 43, 44; averaged", "42, 43, 44; averaged", "42, 43, 44; averaged", "42, 43, 44; averaged", "42, 43, 44; each a separate run"],
    ["Determinism", "reproducible without extra settings", "reproducible (deterministic settings change no value)", "same", "same", "deterministic cuDNN/cuBLAS settings required (see Methods)", "n/a"],
    ["Trainable parameters", *[n_(PC[h]) for h in H[1:]]],
]
s3_md = f"""# Table S3. Architectures, hyperparameters and parameter counts

**Table S3.** Architecture and training settings of TCR-IMAP (pairwise and per-residue variants), the two reimplemented
per-residue comparators (Struct-TCS, Seq-TCS), the reimplemented pairwise comparator (TEIM-Res) and the docking-angle
regressor (AngleNet), as used for every reported (p84v2) result.

{md_table(H, R)}

**Notes.**
- Settings are read from the model and training code; the p84v2 runners execute the same sources with the p84v2 input
  layout (TCR axis 75 rows, segment lengths 7/8/22/6/7/25; MHC axis 84 columns) and change no hyperparameter
  (`results/review/REBUILD_P84V2_REPORT.md` step 3; `CORRECTION_PLAN_P84V2.md` §2).
- Parameter counts were taken by `build_si_tables.py` by instantiating each class on the CPU at the p84v2 layout (no
  forward pass) and counting trainable parameters. The counts for TCR-IMAP (both variants), Struct-TCS and Seq-TCS are
  asserted equal to the values recorded during the build; the TEIM-Res and
  AngleNet counts are recorded here for the first time.
- The pairwise models (TCR-IMAP, TEIM-Res) use no validation subset and no early stopping, so no test complex
  influences their training. The per-residue validation split is described in Methods and evaluated in Figure S4.
- Struct-TCS receives the crystal-structure residue graph of each complex, including test complexes; the other models
  receive amino-acid identities only.
- TEIM-Res: the published model's pretrained peptide autoencoder does not apply to TCR–MHC contacts; all encoders are
  trained from scratch.

**Sources.** `scripts/19_train_pairmap.py` (`PairNet`); `scripts/24b_run_all_p82v1.py` (`fit_predict`: optimizer, loss,
epochs, seeds); `scripts/21_struct_tcs.py` (`OursPerRes`, `TCS`, `GatedGAT`, `AttnBlock`, `GATBlock`);
`scripts/28b_benchmark_p82v1.py` (`fit`, `EPOCHS`, `SEEDS`); `scripts/53_benchmark_p84v2.py`
(`cluster_disjoint_validation`); `scripts/20_build_graphs.py` (graph, `CUT = 8.0`); `scripts/36b_teim_res_p82v1.py`
(`TEIMRes`, training loop); `scripts/53_teim_res_p84v2.py` and `results/review/TEIM_DETERMINISM_REPORT.md`
(determinism); `scripts/16b_docking_angle_p82v1.py` (`AngleNet`, `run_one`); `results/review/VALIDATION_FIX_REPORT.md` §3, §5.
"""
write("S3_hyperparameters", s3_md, {"S3_hyperparameters": pd.DataFrame(R, columns=H)})

# ================================================================= S4 protocols and folds
print("== S4")
proto = json.loads((ROOT / "data/splits_v2/protocols.json").read_text())
A = np.load(T / "tcr_identity_matrix.npy")
pdbs = list(md.pdb); idx = {p: i for i, p in enumerate(pdbs)}
check(A.shape == (174, 174), "S4 identity matrix 174 x 174, metadata order (asserted by the S9 build)")
clu = md.clu70.values
LA = pd.read_csv(T / "leakage_audit_dual.tsv", sep="\t")
PM = pd.read_csv(T / "percomplex_maxsim.tsv", sep="\t")
val = {}
for line in (RV / "validation_fix/B2_split_dryrun.log").read_text().splitlines():
    m = re.match(r"\s*(random|double|allele|temporal|tcr70)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+yes", line)
    if m:
        val[(m[1], int(m[2]))] = dict(test=int(m[3]), train=int(m[4]), val=int(m[5]), val_clusters=int(m[6]), pct=float(m[7]))
check(len(val) == 18, "S4 18 evaluated folds parsed from B2_split_dryrun.log")
s4, s4_csv = [], []
for p in PROTOS:
    for fi, fold in enumerate(proto[p]):
        tei = [idx[x] for x in fold if x in idx]; tri = [i for i in range(174) if i not in set(tei)]
        ms = A[np.ix_(tei, tri)].max(1)
        ctest, ctrain = set(clu[tei]), set(clu[tri])
        la = LA[(LA.protocol == p) & (LA.fold == fi)]
        check(len(la) == 1 and int(la.n_test.iloc[0]) == len(tei) and round(ms.mean(), 3) == float(la.max_sim_mean.iloc[0]),
              f"S4 {p} fold {fi}: n test and mean max similarity = leakage_audit_dual.tsv")
        pmf = PM[(PM.protocol == p) & (PM.fold == fi)]
        if len(pmf):
            check(same(pmf.max_sim.max(), ms.max()), f"S4 {p} fold {fi}: maximum = percomplex_maxsim.tsv")
        ev = len(tei) >= 3
        v = val.get((p, fi))
        if ev:
            check(v is not None and v["test"] == len(tei) and v["train"] + v["val"] == len(tri), f"S4 {p} fold {fi}: per-residue train + validation = non-test complexes (B2 log)")
        row = dict(protocol=p, fold=fi, evaluated=ev, n_test=len(tei), n_train=len(tri), n_test_clusters=len(ctest),
                   n_train_clusters=len(ctrain), n_clusters_shared=len(ctest & ctrain), max_train_test_similarity=float(ms.max()),
                   mean_max_train_test_similarity=float(ms.mean()), frac_test_over_090=float((ms > 0.9).mean()),
                   dup_families_split_allele_dependent=int(la.leak_stored.iloc[0]), dup_families_split_allele_independent=int(la.leak_conservative.iloc[0]),
                   perres_train=(v["train"] if v else np.nan), perres_validation=(v["val"] if v else np.nan),
                   perres_validation_clusters=(v["val_clusters"] if v else np.nan), perres_validation_pct=(v["pct"] if v else np.nan))
        s4_csv.append(row)
        s4.append([plab(p) if fi == 0 else "", fi, "yes" if ev else "**no (n < 3)**", len(tei), len(tri), len(ctest), len(ctrain), len(ctest & ctrain),
                   f"{ms.max():.3f}", f"{ms.mean():.3f}", f"{100 * (ms > 0.9).mean():.1f}", int(la.leak_conservative.iloc[0]), int(la.leak_stored.iloc[0]),
                   (f"{v['train']} / {v['val']} ({v['val_clusters']}; {v['pct']:.1f} %)" if v else "n/a")])
t70 = [r for r in s4_csv if r["protocol"] == "tcr70"]
check(max(r["max_train_test_similarity"] for r in t70) < 0.6958 and round(max(r["max_train_test_similarity"] for r in t70), 3) == 0.696,
      "S4 tcr70 maximum train/test similarity = 0.696 (clustering strings)")
check(all(r["n_clusters_shared"] == 0 for r in t70), "S4 tcr70: no cluster on both sides of any fold")
cov = {p: sum(r["n_test"] for r in s4_csv if r["protocol"] == p and r["evaluated"]) for p in PROTOS}
check(cov == {"random": 174, "double": 30, "allele": 64, "temporal": 60, "tcr70": 174}, f"S4 coverage {cov}")
s4_md = f"""# Table S4. Partitioning protocols: fold composition

**Table S4.** Composition of every fold of the five frozen partitioning protocols: complexes and clu70 clusters on each
side, the train/test TCR similarity, and duplicate families split across the fold boundary.

{md_table(["Protocol", "Fold", "Evaluated", "n test", "n training", "Test clusters", "Training clusters", "Clusters on both sides", "Max train/test similarity", "Mean of per-test maximum", "Test ≥ 0.90 (%)", "Dup. families split (CDR3 pair + peptide)", "Dup. families split (+ allele, stored)", "Per-residue models: training / validation (validation clusters; % held out)"], s4)}

Complexes evaluated: random {cov['random']}, intersection (legacy) {cov['double']}, allele {cov['allele']}, temporal {cov['temporal']},
tcr70 {cov['tcr70']} of 174.

**Notes.**
- n training is every complex outside the test fold (174 − n test). The pairwise models and the frequency priors train on
  all of them. The per-residue models (TCR-IMAP per-residue, Struct-TCS, Seq-TCS) hold out whole clu70 clusters for
  validation until at least 15 % of those complexes are held out; the last column gives that split for evaluated folds.
- Similarity is the Levenshtein similarity ratio (`Levenshtein.ratio`), an edit-distance similarity, not gap-free alignment
  identity, computed on the clustering strings: the six CDR1/CDR2/CDR3 slot segments of both chains, IMGT 27–38, 56–65 and
  105–124 (the last including the J-region framework positions 118–124), one residue per IMGT position. "Max train/test
  similarity" is the largest, over the fold's test complexes, of each test complex's maximum similarity to any training
  complex. Under other string definitions the tcr70 maximum is 0.707 (full CDR-plus-framework sequences; 12 test
  complexes between 0.70 and 0.71, none above 0.80) and 0.682 (canonical CDR sequences alone); see Table S8 and
  `IDENTITY_BOUND_REPORT.md`.
- **Allele and intersection coverage limitation.** The allele and intersection protocols cover only 64 and 30 of the 174 complexes.
  Allele annotation is regex-parsed from the free-text `compound` field of the STCRDab summary, so 108 complexes have
  none. In `scripts/23_make_protocols.py` (lines 26 and 33) the allele groups for the folds were formed with
  `fillna("NA")`, but fold membership was then tested with `.isin()` on the original column, where those 108 complexes
  are still missing; they therefore enter no allele or intersection test fold. Allele fold 3 holds only two complexes and is
  dropped by the three-complex minimum used in every training script. The frozen protocol file was not regenerated.
- Duplicate families: identical CDR3α–CDR3β pair and peptide (allele-independent; the definition reported) and, as
  stored in `complex_metadata.tsv`, identical CDR3 pair, allele and peptide (allele-dependent; grouped over 108 missing
  alleles, so it both over- and under-merges). tcr70 is 0 in every fold under both.
- Fold-level duplicate counts include allele fold 3; the mean over the four evaluated allele folds is 1.25 (Figure 2E),
  and 1.2 over all five.

**Sources.** `data/splits_v2/protocols.json` (SHA-256 verified), `data/processed/complex_metadata.tsv` (clu70; SHA-256
verified), `results/tables/tcr_identity_matrix.npy` (similarity), checked against `results/tables/leakage_audit_dual.tsv`
(n test, mean of per-test maximum; duplicate families taken from its `leak_conservative` and `leak_stored` columns) and
`results/tables/percomplex_maxsim.tsv` (per-fold maximum); per-residue validation split from
`results/review/validation_fix/B2_split_dryrun.log`. Computed by `build_si_tables.py` (arithmetic only).
"""
write("S4_protocol_folds", s4_md, {"S4_protocol_folds": pd.DataFrame(s4_csv)})

# ================================================================= S5 sweep
print("== S5")
_frozen_order = pd.read_csv(T / "similarity_curve_p84v2.csv")                   # 2,000, seed 0 (superseded)
sc = (pd.read_csv(T / "similarity_curve_p84v2_authoritative.csv")                # 20,000 resamples, seed 42
        .set_index("threshold").loc[_frozen_order.threshold].reset_index())      # keep the frozen row order
sc_old = _frozen_order.set_index("threshold")
check(bool(np.allclose(sc["diff"].values, _frozen_order["diff"].values)),
      "S5 point estimates identical between the authoritative and superseded bootstrap settings")
comp = pd.read_csv(T / "similarity_curve_composition_p84v2.csv")
cf = pd.read_csv(T / "similarity_curve_composition_folds_p84v2.csv")
check(list(sc.threshold) == list(comp.threshold), "S5 thresholds align")
s5, s5_csv = [], []
for r, c in zip(sc.itertuples(), comp.itertuples()):
    folds = cf[cf.threshold == r.threshold].sort_values("fold")
    check(int(folds.n_train_kept.min()) == c.n_train_min and int(folds.n_train_kept.max()) == c.n_train_max and int(c.n_test) == r.n,
          f"S5 threshold {r.threshold}: fold table = summary table, n test = sweep n")
    thr = "none" if r.threshold > 1 else f"< {r.threshold:.2f}"
    s5.append([thr, r.n, " / ".join(str(int(x)) for x in folds.n_train_kept), f"{c.n_train_mean:.1f}", f"{c.achieved_max_sim_max:.3f}",
               f"{c.achieved_max_sim_mean:.3f}", f"{r.mean_max_sim:.4f}", fx(r.model), fx(r.freq), sg(r.diff),
               ci(r.lo, r.hi), pstr(r.p),
               ci(sc_old.loc[r.threshold].lo, sc_old.loc[r.threshold].hi), pstr(sc_old.loc[r.threshold].p)])
    s5_csv.append(dict(threshold=thr, n_test=r.n, n_train_per_fold=" / ".join(str(int(x)) for x in folds.n_train_kept),
                       n_train_mean=c.n_train_mean, achieved_max_similarity=c.achieved_max_sim_max, achieved_mean_similarity=c.achieved_max_sim_mean,
                       mean_max_sim_column_upper_bound=r.mean_max_sim, tcrimap_macro_auprc=r.model, prior_macro_auprc=r.freq,
                       diff=r.diff, ci_lo=r.lo, ci_hi=r.hi, p=r.p, resamples="20,000 (seed 42)",
                       superseded_ci_lo=sc_old.loc[r.threshold].lo, superseded_ci_hi=sc_old.loc[r.threshold].hi,
                       superseded_p=sc_old.loc[r.threshold].p, superseded_resamples="2,000 (seed 0)"))
s5_md = f"""# Table S5. Identity-threshold sweep: composition and results

**Table S5.** For each threshold on the TCR Levenshtein similarity ratio, the training complexes kept in each of the
five folds, the train/test similarity actually achieved, and pairwise macro-AUPRC of TCR-IMAP and the pair-frequency
prior with their paired difference.

{md_table(["Threshold", "n test", "n training per fold (folds 0–4)", "Mean n training", "Achieved max similarity", "Achieved mean similarity", "`mean_max_sim` column (upper bound)", "TCR-IMAP", "Prior", "Δ", "95 % CI", "p", "95 % CI (superseded)", "p (superseded)"], s5)}

**Notes.**
- The sweep uses one fixed random five-fold partition (`np.random.default_rng(0)`, independent of the frozen
  protocols; folds of 35, 35, 35, 35 and 34). At each threshold every non-test complex whose similarity to **any** test
  complex of its fold reaches the threshold is removed from that fold's training set; a fold is used only if at least 30
  training complexes remain (all five are used at every threshold). The test composition never changes: all 174
  complexes are evaluated at every threshold. The model and prior are retrained at every threshold.
- Achieved similarity: for each test complex, its maximum similarity to the training complexes its fold kept; the
  maximum and the mean over the 174 test complexes.
- **The `mean_max_sim` column recorded with the sweep is not the achieved similarity.** For each test complex it takes
  the maximum similarity to any *other complex of the 174* whose similarity is below the threshold, including complexes
  in the same test fold and training complexes the filter removed, and averages over test complexes. It is therefore
  an upper bound on the achieved mean (here 0.008–0.029 above it, and different even with no restriction, because it
  also counts same-fold neighbors).
- Δ = TCR-IMAP − prior. The reported **95 % CI and p use the project standard: paired clu70 cluster bootstrap,
  20,000 resamples, seed 42**, p floored at 1/20,000 and reported as < 0.001 at the floor.
- The two rightmost columns carry the **superseded** values first recorded for this stage, which used the
  2,000-resample seed-0 default. The sweep originally stored no per-complex arrays, so it could not be recomputed at
  the standard; those arrays were later regenerated under the stage's own configuration and seeds, each threshold
  reproducing the frozen aggregate to four decimal places, which is what allows the recomputation.
  **The point estimates are identical under both settings** (asserted in the build). The intervals and p values differ
  only by Monte Carlo error in the bootstrap, which 2,000 resamples do not resolve near an interval endpoint; at
  < 0.70 that moves the interval endpoint across zero (p 0.052 → 0.041). This is a numerical property of the
  resampling, not a second analysis and not evidence about evaluation design.
- Similarity uses the clustering strings (see Table S4).

**Sources.** `results/tables/similarity_curve_p84v2.csv` (results, `mean_max_sim`),
`results/tables/similarity_curve_composition_p84v2.csv` and `similarity_curve_composition_folds_p84v2.csv` (composition;
built by `scripts/57_sweep_composition_p84v2.py` from `tcr_identity_matrix.npy` and the sweep's own fold seed;
`results/review/SI_FIGURES_REPORT.md`, section S9). Sweep procedure: `scripts/25b_similarity_curve_p82v1.py`, executed on
the p84v2 inputs by `scripts/53_similarity_curve_p84v2.py`.
"""
write("S5_identity_sweep", s5_md, {"S5_identity_sweep": pd.DataFrame(s5_csv)})

# ================================================================= S6 MD selection and QC
print("== S6")
SYS = ["5yxn", "7n1f", "3vxs", "8eo8"]
# Selection record, md/MD_SYSTEM_SELECTION_CRITERIA.md section 8 (p82v1 values, as committed before any simulation)
SEL = {"5yxn": ("high-Δ", 0.2503, 0.5360, 0.2856, 2.03, 3, "HLA-A*02:01", "KLVALGINAV", "pre-identified anchor"),
       "3vxs": ("high-Δ", 0.2167, 0.6196, 0.4029, 1.80, 28, "HLA-A*24:01", "RYPLTLGWCF", "next-ranked eligible"),
       "7n1f": ("low-Δ", -0.2314, 0.4490, 0.6804, 2.39, 42, "not annotated", "YLQPRTFLL", "pre-identified anchor"),
       "8eo8": ("low-Δ", -0.1946, 0.1695, 0.3641, 2.30, 38, "HLA-B*35:01", "LPFDKATIM", "next-ranked eligible")}
RANK_REC = {"5yxn": (12, 0.245), "3vxs": (4, 0.266), "7n1f": (167, -0.111), "8eo8": (154, -0.065)}
ptr = pt.copy(); ptr["rank"] = ptr.delta.rank(ascending=False, method="first").astype(int)
check(np.allclose(ptr.delta, ptr.grip - ptr.prior), "S6 percomplex_tcr70_paired delta = grip − prior")
mdi = md.set_index("pdb")
s6a, s6a_csv = [], []
for s in ["5yxn", "3vxs", "7n1f", "8eo8"]:
    st, d82, m82, p82, res_, c70, al, pep, how = SEL[s]
    r = ptr[ptr.pdb == s].iloc[0]
    # tolerance 0.001: the revision records 7n1f as -0.111, a double rounding of -0.110461 (-> -0.1105 -> -0.111)
    check(int(r["rank"]) == RANK_REC[s][0] and abs(float(r.delta) - RANK_REC[s][1]) < 0.001,
          f"S6 {s} p84v2 rank {int(r['rank'])} and Δ {float(r.delta):+.6f} = criteria document revision of 24 Oct (±0.001)")
    check(int(mdi.loc[s, "clu70"]) == c70 and abs(float(mdi.loc[s, "resolution"]) - res_) < 0.006, f"S6 {s} clu70 and resolution = complex_metadata")
    s6a.append([s, st, how, sg(d82), fx(m82), fx(p82), f"{res_:.2f}", c70, al, pep, f"{int(r['rank'])} / 174", sg(float(r.delta)), fx(float(r.grip)), fx(float(r.prior))])
    s6a_csv.append(dict(pdb=s, stratum=st, selected_as=how, p82v1_delta=d82, p82v1_tcrimap=m82, p82v1_prior=p82, resolution_A=res_, clu70=c70,
                        allele=al, peptide=pep, p84v2_rank_of_174=int(r["rank"]), p84v2_delta=float(r.delta), p84v2_tcrimap=float(r.grip), p84v2_prior=float(r.prior)))

qc = {}
cur = None
for line in (RV / "characterisation_fix/S2_qc_all12.log").read_text().splitlines():
    m = re.match(r"#+ (\w{4}) rep(\d) #+", line)
    if m:
        cur = (m[1], int(m[2])); qc[cur] = {"if": {}}; continue
    if cur is None:
        continue
    q = qc[cur]
    if m := re.search(r"checkpoint step : (\d+) / (\d+)", line): q["step"] = int(m[1])
    if m := re.search(r"queue events\s+: (.+?)\s+\| mdrun sessions (\d+) \| LINCS warnings (\d+)", line):
        q["events"], q["sessions"], q["lincs"] = m[1].strip(), int(m[2]), int(m[3])
    if m := re.search(r"trajectory\s+: (\d+) frames every (\d+) ps.*irregular steps (\d+)", line): q["frames"], q["irr"] = int(m[1]), int(m[3])
    for key in ("Temperature", "Pressure", "Density"):
        if m := re.search(rf"{key}\s+: avg\s+([-\d.]+)\s+sd\s+([-\d.]+)\s+drift\(20-200 ns\)\s+([-+\d.]+)", line):
            q[key] = (float(m[1]), float(m[2]), float(m[3]))
    if m := re.search(r"backbone RMSD\s+: mean\(100-200 ns\) ([\d.]+) nm, max ([\d.]+), final ([\d.]+)", line):
        q["qc_rmsd"] = tuple(float(x) for x in m.groups())
    if m := re.match(r"\s+(pep-MHC|pep-TCRa|pep-TCRb|TCRa-MHC|TCRb-MHC|MHC-b2m|TCRa-TCRb)\s+(\d+) ->\s+(\d+) ->\s+(\d+)", line):
        q["if"][m[1]] = tuple(int(x) for x in m.groups()[1:])
check(len(qc) == 12, "S6 12 replicates parsed from S2_qc_all12.log")
RET_REC = {("3vxs", 1): 71, ("3vxs", 2): 85, ("3vxs", 3): 82, ("7n1f", 1): 58, ("7n1f", 2): 57, ("7n1f", 3): 74,
           ("8eo8", 1): 73, ("8eo8", 2): 58, ("8eo8", 3): 75, ("5yxn", 1): 21, ("5yxn", 2): 29, ("5yxn", 3): 34}
cs = pd.read_csv(T / "md_characterisation_summary.csv", float_precision="round_trip")
ts = pd.read_csv(T / "md_characterisation_timeseries.csv", float_precision="round_trip")
dfr = pd.read_csv(T / "md_occupancy_docked_frames.csv", float_precision="round_trip")
s6b, s6b_csv, RET_DIFF = [], [], []
for s in SYS:
    for rep in (1, 2, 3):
        q = qc[(s, rep)]
        start = sum(q["if"][k][0] for k in ("pep-TCRa", "pep-TCRb", "TCRa-MHC", "TCRb-MHC"))
        end = sum(q["if"][k][2] for k in ("pep-TCRa", "pep-TCRb", "TCRa-MHC", "TCRb-MHC"))
        ret = 100 * end / start
        check(abs(ret - RET_REC[(s, rep)]) <= 1.5, f"S6 {s} rep{rep} TCR–pMHC retention {ret:.1f} % = criteria document table ({RET_REC[(s, rep)]} %) within 1.5 points")
        if round(ret) != RET_REC[(s, rep)]:
            RET_DIFF.append(f"{s} rep{rep} {ret:.1f} % (recorded {RET_REC[(s, rep)]} %)")
        c = cs[(cs.pdb == s) & (cs.rep == rep)].iloc[0]
        # the timeseries also holds 181 SASA-only rows per replicate (1 ns stride, backbone NaN); keep the RMSD frames
        tt = ts[(ts.pdb == s) & (ts.rep == rep) & (ts.time_ps >= 20000) & (ts.time_ps <= 200000) & ts.backbone_rmsd.notna()]
        check(len(tt) == 9001 and same(tt.backbone_rmsd.mean(), c.backbone_rmsd_mean, 1e-9), f"S6 {s} rep{rep}: 9,001 window frames; mean backbone RMSD = summary table")
        dk = dfr[(dfr.pdb == s) & (dfr.rep == rep)].iloc[0]
        check(q["step"] == 50_000_000 and q["frames"] == 10001 and q["irr"] == 0 and q["lincs"] == 0, f"S6 {s} rep{rep}: complete, 10,001 frames, 0 irregular, 0 LINCS")
        Tm, P, D = q["Temperature"], q["Pressure"], q["Density"]
        s6b.append([s, rep, n_(q["step"]), n_(q["frames"]), q["lincs"], q["sessions"], f"{Tm[0]:.2f} ± {Tm[1]:.2f}", fx(P[0], 2) + f" ± {P[1]:.1f}",
                    f"{D[0]:.2f} ({sg(D[2], 2)})", f"{c.backbone_rmsd_mean:.2f}", f"{tt.backbone_rmsd.max():.2f}", f"{start} → {end} ({ret:.1f} %)",
                    f"{100 * dk.docked_fraction_count:.0f} / {100 * dk.docked_fraction_identity:.0f}"])
        s6b_csv.append(dict(pdb=s, rep=rep, checkpoint_step=q["step"], frames_20ps=q["frames"], irregular_steps=q["irr"], lincs_warnings=q["lincs"],
                            mdrun_sessions=q["sessions"], queue_events=q["events"], temperature_K_mean=Tm[0], temperature_K_sd=Tm[1],
                            pressure_bar_mean=P[0], pressure_bar_sd=P[1], density_kgm3_mean=D[0], density_drift_20_200ns=D[2],
                            backbone_rmsd_A_mean_20_200ns=float(c.backbone_rmsd_mean), backbone_rmsd_A_max_20_200ns=float(tt.backbone_rmsd.max()),
                            qc_backbone_rmsd_nm_mean_100_200ns=q["qc_rmsd"][0], qc_backbone_rmsd_nm_max=q["qc_rmsd"][1],
                            tcr_pmhc_contacts_t0=start, tcr_pmhc_contacts_200ns=end, tcr_pmhc_retention_pct=ret,
                            docked_pct_count=100 * dk.docked_fraction_count, docked_pct_identity=100 * dk.docked_fraction_identity))
rm = pd.DataFrame(s6b_csv).groupby("pdb").backbone_rmsd_A_mean_20_200ns.mean()
check(all(round(float(rm[s]), 2) == v for s, v in (("5yxn", 4.25), ("7n1f", 3.98), ("3vxs", 3.12), ("8eo8", 3.80))), "S6 system mean backbone RMSD = Figure 7E caption")
s6_md = f"""# Table S6. Molecular dynamics: system selection and per-replicate quality control

**Table S6.** (A) The four simulated TCR–pMHC complexes: the predeclared selection basis, structure and their standing
under the corrected (p84v2) analysis. (B) Quality control of each of the twelve 200 ns production replicates.

## A. Systems

{md_table(["PDB", "Stratum", "Selected as", "Δ (p82v1)", "TCR-IMAP (p82v1)", "Prior (p82v1)", "Resolution (Å)", "clu70", "Allele", "Peptide", "p84v2 rank", "Δ (p84v2)", "TCR-IMAP (p84v2)", "Prior (p84v2)"], s6a)}

## B. Production replicates

{md_table(["PDB", "Rep", "Steps", "Frames (20 ps)", "LINCS warnings", "mdrun sessions", "T (K), mean ± SD", "P (bar), mean ± SD", "Density (kg m⁻³), mean (drift 20–200 ns)", "Backbone RMSD mean (Å)", "Backbone RMSD max (Å)", "TCR–pMHC contacts t = 0 → 200 ns (retained)", "Docked % (count / identity)"], s6b)}

**Notes.**
- Δ = per-complex pairwise AUPRC of TCR-IMAP − pair-frequency prior under tcr70. **Selection used the p82v1 values**
  (the analysis version current when the criteria were committed), mechanically, before any system was built or
  simulated: the two highest- and two lowest-Δ complexes passing the predeclared filters (resolution ≤ 2.5 Å; distinct
  clu70 clusters; complete α1/α2 window; a single peptide of standard residues, no covalent TCR–peptide link, engineered
  interchain disulfide or single-chain construct; member of the 174-complex set). 5yxn and 7n1f were pre-identified
  anchors, included because they passed every filter. The selection was not revisited after the p84v2 correction; the
  p84v2 rank (1 = largest Δ) is given for transparency. The p84v2 extremes are 3qfj (+0.279) and 3qeq (−0.158).
  7n1f's p84v2 Δ is −0.110 (stored −0.110461). The criteria document's revision first gave −0.111, a double rounding
  through −0.1105, and was corrected to −0.110 on 28 Oct 2026; the rank, 167, is unchanged.
- Every replicate: 50,000,000 steps of 4 fs (200 ns), checkpoint complete, 10,001 frames at 20 ps with no irregular
  timesteps, 0 LINCS warnings. Interruptions were resumed from checkpoint, which the protocol allows. Temperature,
  pressure and density are the QC script's averages over the production run, with the density drift over 20–200 ns.
- Backbone RMSD: protein Cα RMSD from the energy-minimized structure over the analyzed window (20–200 ns, every frame,
  9,001 frames), on trajectories made whole and jump-free (the characterization recomputed after the periodic-boundary
  fix). The QC script's own coarser RMSD (every 10 ns, nm) is kept in the CSV.
- TCR–pMHC contacts: heavy-atom pairs < 4.5 Å between TCRα/TCRβ and the peptide plus MHC heavy chain (β2m excluded),
  summed over the four interfaces, at t = 0 and 200 ns; retained = 200 ns / t = 0, computed here from the QC log and
  given to one decimal. The criteria document's retention table gives the same values: its replicate-2 entries,
  recorded as integers in QC round 2 (17 Oct), were corrected to these log values on 28 Oct 2026.
- Docked %: share of the 9,001 window frames counted as docked, under the count reading (at least half as many TCR–pMHC
  heavy-atom contacts as at t = 0; the predeclared rule) and the identity reading (at least half of the specific t = 0
  contact pairs present).
- mdrun sessions are as logged by the QC script. They differ from the narrative in the criteria document ("All 12
  replicates against the predeclared protocol"), which names interruptions only for 5yxn rep1, 5yxn rep3, 7n1f rep2 and
  3vxs rep2; the log also records 3 sessions for 5yxn rep2 and 2 for 7n1f rep1. No other quantity is affected.

**Sources.** A: `md/MD_SYSTEM_SELECTION_CRITERIA.md` §8 (selection record: p82v1 Δ, TCR-IMAP, prior, resolution, clu70,
allele, peptide) and the revision in which the analysis moves to the p84v2 pair universe (p84v2 ranks); p84v2 rank and values
recomputed from `results/tables/percomplex_tcr70_paired_p84v2.tsv` and asserted equal to that revision; clu70 and
resolution checked against `data/processed/complex_metadata.tsv`. B: `results/review/characterisation_fix/S2_qc_all12.log`
(`md/qc_replicate.sh` on all 12 replicates: steps, frames, LINCS, sessions, thermodynamics, interface contacts; retention
checked against the criteria document's retention table), `results/tables/md_characterisation_summary.csv` and
`md_characterisation_timeseries.csv` (backbone RMSD), `results/tables/md_occupancy_docked_frames.csv` (docked %).
"""
write("S6_md_selection_qc", s6_md, {"S6_md_systems": pd.DataFrame(s6a_csv), "S6_md_replicates": pd.DataFrame(s6b_csv)})

# ================================================================= S7 MD occupancy statistics
print("== S7")
rd = lambda f: pd.read_csv(T / f, float_precision="round_trip")  # noqa: E731
rho, strata, topk, ck1, dr, ss = (rd(f) for f in ("md_occupancy_rho.csv", "md_occupancy_strata.csv", "md_occupancy_topk.csv",
                                                   "md_occupancy_check1_positive.csv", "md_occupancy_docked_rho.csv", "md_score_score_spearman.csv"))
SC = ["TCR-IMAP", "pair-frequency prior"]


def ptxt(p):
    return "n/a" if pd.isna(p) else ("< 0.001" if p <= 0.0005 else f"{p:.3f}")


s7a = []
for lvl, pdb_ in [("pooled", "all four")] + [("system", s) for s in SYS]:
    a, b = (rho[(rho.score == sc_) & (rho.level == lvl) & (rho.pdb == pdb_)].iloc[0] for sc_ in SC)
    check(a.n_crystal_negative == b.n_crystal_negative and a.n_nonzero_occupancy == b.n_nonzero_occupancy, f"S7 {pdb_}: both scores on the same pairs")
    s7a.append([pdb_ if lvl == "system" else "**pooled (all four)**", "all", n_(a.n_crystal_negative), n_(a.n_nonzero_occupancy), n_(a.n_occupancy_ge_025),
                f"{a.rho:.3f}", ptxt(a.perm_p_approximate), f"{b.rho:.3f}", ptxt(b.perm_p_approximate)])
    for stv in ("4.5-6 A", "6-8 A", "> 8 A"):
        a, b = (strata[(strata.score == sc_) & (strata.level == lvl) & (strata.pdb == pdb_) & (strata.stratum == stv)].iloc[0] for sc_ in SC)
        s7a.append(["", stv.replace("A", "Å").replace("-", "–"), n_(a.n), n_(a.n_nonzero_occupancy), n_(a.n_occupancy_ge_025),
                    f"{a.rho:.3f}", ptxt(a.perm_p_approximate), f"{b.rho:.3f}", ptxt(b.perm_p_approximate)])
check(bool(rho.interpretable.all()) and bool(strata.interpretable.all()), "S7 every primary and stratum row is interpretable")
s7b = []
for s in SYS:
    for rep in (1, 2, 3):
        a, b = (rho[(rho.score == sc_) & (rho.level == "replicate") & (rho.pdb == s) & (rho.rep.astype(str) == str(rep))].iloc[0] for sc_ in SC)
        s7b.append([s, rep, n_(a.n_crystal_negative), n_(a.n_nonzero_occupancy), n_(a.n_occupancy_ge_025), f"{a.rho:.3f}", f"{b.rho:.3f}"])
s7c = []
for s in SYS:
    for sc_ in SC:
        r = topk[(topk.score == sc_) & (topk.pdb == s)].iloc[0]
        s7c.append([s if sc_ == SC[0] else "", sc_, r.k, f"{int(r.n_top_k_ge)} / {r.k} ({r.frac_top_k:.3f})", f"{r.base_rate:.4f}",
                    f"{r.frac_top_k / r.base_rate:.0f}×", "yes" if r.tie_at_k else "no"])
enr = topk.frac_top_k / topk.base_rate
check(round(enr.min()) == 28 and round(enr.max()) == 88 and topk.n_top_k_ge.min() == 1 and topk.n_top_k_ge.max() == 7, "S7 top-k enrichment 28–88× on 1–7 pairs (recorded value)")
s7d = [[r.pdb, r.rep, r.n, f"{r['mean']:.3f}", f"{r['median']:.3f}", f"{r['min']:.3f}", f"{r.n_ge_meaningful} / {r.n}"] for _, r in ck1.iterrows()]
s7e = []
for reading in ("count", "identity"):
    for _, a in dr[(dr.score == SC[0]) & (dr.reading == reading)].iterrows():
        b = dr[(dr.score == SC[1]) & (dr.reading == reading) & (dr.level == a.level) & (dr.pdb == a.pdb)].iloc[0]
        nc = "n/a" if pd.isna(a.n_replicates_contributing) else int(a.n_replicates_contributing)
        s7e.append([reading, a.pdb, nc, n_(a.n_crystal_negative), "n/a" if pd.isna(a.n_nonzero_occupancy) else n_(a.n_nonzero_occupancy),
                    "n/a" if pd.isna(a.rho) else f"{a.rho:.3f}", "n/a" if pd.isna(b.rho) else f"{b.rho:.3f}", "yes" if a.interpretable else "no (no contributing replicate)"])
s7f = [[r.pdb, n_(r.n_crystal_negative), f"{r.spearman_rho:.3f}"] for _, r in ss.iterrows()]
s7_md = f"""# Table S7. Molecular dynamics: occupancy statistics

**Table S7.** Spearman correlation between score and solution-phase contact occupancy among crystal-negative residue
pairs, for TCR-IMAP and the pair-frequency prior (the predeclared control), with the occupancy counts each correlation
rests on. (A) Pooled and per system, overall and by crystal-frame distance stratum. (B) Per replicate. (C) Top-k
enrichment. (D) Occupancy of crystal-positive pairs (check 1). (E) Docked-frame secondary analysis. (F) Correlation
between the two scores.

## A. Primary statistic, pooled and per system, by stratum

{md_table(["System", "Stratum", "n pairs", "n non-zero", "n ≥ 0.25", "ρ TCR-IMAP", "p (approx.)", "ρ prior", "p (approx.)"], s7a)}

## B. Per replicate (all crystal-negative pairs)

{md_table(["System", "Rep", "n pairs", "n non-zero", "n ≥ 0.25", "ρ TCR-IMAP", "ρ prior"], s7b)}

## C. Top-k enrichment (k = the system's number of crystal-positive pairs)

{md_table(["System", "Score", "k", "Top-k pairs with occupancy ≥ 0.25 (fraction)", "Base rate (all crystal-negative pairs)", "Enrichment", "Tie in score at rank k"], s7c)}

## D. Crystal-positive occupancy (check 1)

{md_table(["System", "Rep", "n crystal-positive pairs", "Mean", "Median", "Min", "Occupancy ≥ 0.25"], s7d)}

## E. Docked-frame secondary analysis (system occupancy over docked frames of contributing replicates)

{md_table(["Reading", "System", "Replicates contributing", "n pairs", "n non-zero", "ρ TCR-IMAP", "ρ prior", "Interpretable"], s7e)}

## F. Score–score Spearman correlation among crystal-negative pairs

{md_table(["System", "n pairs", "ρ (TCR-IMAP vs prior)"], s7f)}

**Notes.**
- Occupancy of a pair: fraction of the 9,001 frames between 20 and 200 ns in which its minimum heavy-atom distance is
  ≤ 4.5 Å (the crystal-label definition); system occupancy is the mean of the three replicates. Pair universe: resolved
  pairs of the p84v2 label array (75 TCR rows × 84 MHC columns); crystal-negative = not a contact in the crystal
  structure. Scores: TCR-IMAP's out-of-fold tcr70 pairwise predictions and the within-fold pair-frequency prior.
- p: within-system permutation test, 10,000 permutations of the scores among a system's crystal-negative pairs (within
  system for the pooled rows), seed 42, two-sided, **approximate** because pairs within a system are not independent;
  p is floored at 1/10,001 and shown as < 0.001. Per-replicate rows carry no p. No test of the difference between the
  two scores was made.
- Strata: the pair's minimum heavy-atom distance in the crystal structure, 4.5–6, 6–8 and > 8 Å (fixed before the run).
- Every row in A and B rests on at least 10 non-zero occupancies and is flagged `interpretable`; in E, the identity
  reading has no contributing replicate for 5yxn and 8eo8 (no rho; reported, not interpreted).
- C: enrichment = fraction of top-k pairs with occupancy ≥ 0.25 ÷ base rate; it rests on 1–7 pairs per system. A tie in
  score at rank k is flagged; the top-k set uses a stable descending sort in row-major pair order.
- D: no numeric acceptance criterion was predeclared for check 1.
- E: docked frame, count reading (predeclared, primary for this analysis): at least half as many TCR–pMHC heavy-atom
  contacts as at t = 0; identity reading: at least half of the specific t = 0 contact pairs present. A replicate
  contributes only with ≥ 1,000 docked frames.

**Sources.** `results/tables/md_occupancy_rho.csv` (A overall rows, B), `md_occupancy_strata.csv` (A strata),
`md_occupancy_topk.csv` (C), `md_occupancy_check1_positive.csv` (D), `md_occupancy_docked_rho.csv` (E),
`md_score_score_spearman.csv` (F); all written once by the single authorised run (manifest `md_occupancy_manifest.csv`)
and read with `float_precision="round_trip"`. Definitions: `md/MD_SYSTEM_SELECTION_CRITERIA.md` (analysis definitions of
20 Oct; revisions and decisions of 24 Oct).
"""
write("S7_md_occupancy", s7_md, {"S7_md_rho": rho, "S7_md_strata": strata, "S7_md_topk": topk, "S7_md_check1": ck1,
                                 "S7_md_docked_rho": dr, "S7_md_score_score": ss})

# ================================================================= S8 dataset descriptors
print("== S8")
D = A.copy(); np.fill_diagonal(D, -1); nnb = D.max(1)
twins = int((nnb >= 0.999).sum())
pairs = [(i, j) for i in range(174) for j in range(i + 1, 174) if A[i, j] >= 0.999]
dp = sum(md.peptide[i] != md.peptide[j] for i, j in pairs)
cnt = md.dup_family.value_counts(); keep = md.dup_family.map(cnt).values == 1
Ak = A[np.ix_(keep, keep)].copy(); np.fill_diagonal(Ak, -1); med = float(np.median(Ak.max(1)))
H6 = pd.read_csv(RV / "identity_bound/handoff6_statements.tsv", sep="\t")
h = H6[H6.definition.str.startswith("(1)")].set_index("statement").value
check(h["complexes with an exact TCR twin (NN >= 0.999)"] == f"{twins}/174 = {100 * twins / 174:.1f} %", f"S8 exact twins {twins}/174 = handoff6_statements.tsv")
check(h["twin pairs differing in peptide"] == f"{dp}/{len(pairs)} = {100 * dp / len(pairs):.1f} %", f"S8 twin pairs differing in peptide {dp}/{len(pairs)} = handoff6_statements.tsv")
check(h["deduplication: complexes in multi-member (allele-dependent) families removed -> remaining, median NN"] ==
      f"{int((~keep).sum())} removed -> {int(keep.sum())}, median NN {med:.3f}", f"S8 deduplication {int((~keep).sum())} removed -> {int(keep.sum())}, median NN {med:.3f} = handoff6_statements.tsv")
LS = pd.read_csv(T / "gdomain182_region_loop_summary_p84v2.csv")
helix = int(LS[LS.region.isin(["alpha1 helix", "alpha2 helix"])].n_contacts.sum()); tot = int(LS.n_contacts.sum())
floor_ = LS[LS.region == "beta-sheet floor"]
check((helix, tot) == (4674, 5294), f"S8 helix contacts {helix} of {tot} = REBUILD_P84V2_REPORT.md step 4 (4,674 / 5,294)")
check(int(floor_.n_contacts.sum()) == 0 and int(floor_.n_pairs_observed.sum()) == 595952, "S8 β-sheet floor: 0 contacts in 595,952 resolved pairs")
s8 = [
    ["Exact TCR twins", f"{twins} of 174 complexes ({100 * twins / 174:.1f} %)",
     "Complexes whose nearest neighbor (another complex) has a Levenshtein similarity ratio ≥ 0.999 on the clustering strings.",
     "`tcr_identity_matrix.npy`; `identity_bound/handoff6_statements.tsv` (script `54b_dupfamily_and_handoff6.py`)"],
    ["Twin pairs differing in peptide", f"{dp} of {len(pairs)} twin pairs ({100 * dp / len(pairs):.1f} %)",
     "Among all complex pairs with similarity ≥ 0.999, those whose `peptide` field differs.",
     "`tcr_identity_matrix.npy`, `complex_metadata.tsv`; `handoff6_statements.tsv` (54b)"],
    ["Deduplication is insufficient", f"{int((~keep).sum())} complexes removed → {int(keep.sum())} remain; median nearest-neighbor similarity {med:.3f}",
     "Remove every complex in a duplicate family with more than one member (the stored `dup_family`: identical CDR3 pair, allele and peptide); for each remaining complex, its maximum similarity to any other remaining complex; median over complexes.",
     "`complex_metadata.tsv` (`dup_family`), `tcr_identity_matrix.npy`; `handoff6_statements.tsv` (54b)"],
    ["Helix share of G-domain contacts", f"{helix:,} of {tot:,} contacts ({100 * helix / tot:.1f} %)",
     "TCR–MHC contacts (minimum heavy-atom distance ≤ 4.5 Å) over the 182-position α1/α2 G-domain, summed over the six loops, on positions whose DSSP majority class is α1 or α2 helix.",
     "`gdomain182_region_loop_summary_p84v2.csv` (script `53_gdomain182_descriptive_p84v2.py`); `REBUILD_P84V2_REPORT.md` step 4"],
    ["β-sheet floor", f"{int(floor_.n_contacts.sum())} contacts in {int(floor_.n_pairs_observed.sum()):,} resolved pairs",
     "Same table, positions classified β-sheet floor (56 positions).", "`gdomain182_region_loop_summary_p84v2.csv`"],
]
s8_md = f"""# Table S8. Dataset descriptors

**Table S8.** Descriptive statistics of the 174-complex benchmark that are quoted in the text but were previously
recorded only in review documents, each with the computation that produces it.

{md_table(["Quantity", "Value", "Computation", "Source"], s8)}

**Notes.**
- Similarity: Levenshtein similarity ratio on the clustering strings (CDR1/CDR2/CDR3 slot segments of both chains, IMGT
  27–38, 56–65 and 105–124, FR4 118–124 included, one residue per IMGT position), the strings that define clu70 and the
  frozen splits. Under the other string definitions of `IDENTITY_BOUND_REPORT.md` the twin and peptide statistics are
  unchanged and the deduplication median is 0.968 (full CDR-plus-framework sequences) or 0.958 (canonical CDRs alone).
- The deduplication uses the stored, allele-dependent `dup_family` (CDR3 pair, allele, peptide), as the original
  statement did; its grouping over 108 missing alleles is a known limitation of the frozen metadata (Table S4).
- The helix share is computed over the 182-position G-domain (IMGT 1–90 and 1001–1090 plus 1061A and 1072A), of whose
  5,294 contacts the 84-column window holds 5,281. Region classes are DSSP majority calls (helix H/G/I, sheet E/B,
  ≥ 30 resolving structures, majority > 0.5).
- Every value was recomputed by `build_si_tables.py` from the committed inputs and asserted equal to the recorded value.

**Sources.** Listed per row. Recorded values: `results/review/identity_bound/handoff6_statements.tsv` (definition 1),
`results/review/REBUILD_P84V2_REPORT.md` step 4; the remaining recorded values are those asserted during the build.
"""
write("S8_dataset_descriptors", s8_md, {"S8_dataset_descriptors": pd.DataFrame(s8, columns=["quantity", "value", "computation", "source"])})

# ================================================================= INDEX
caps = {
    "S1_benchmark": ("Full benchmark: every method under every protocol",
                     "Per-complex macro-AUPRC and ROC AUC of every method under every protocol, n complexes and clusters, paired differences from TCR-IMAP (per-residue, pairwise, comparator subset, AngleNet)",
                     "`percomplex_{protocol}_p84v2.npz`, `benchmark_p84v2.csv`, `protocol_main_p84v2.csv`, `pairwise_ours_p84v2.npz`, `teim_res_pairwise_p84v2.npz`, `structcs_setting_p84v2_*.csv`, `angle_runs_p84v2.jsonl`; `figures_p84v2/stats_20k_seed42.tsv`, `validation_fix/F_structcs_differences_20k.tsv`"),
    "S2_ablations": ("Ablations of TCR-IMAP (tcr70, pairwise task)",
                     "All eight ablations: Δ, 95 % CI, Bonferroni-simultaneous 99.375 % interval, raw and Holm-adjusted p, survival",
                     "`_partial/p84v2/29b_ablations/*.npz`, `ablations_p84v2.csv`; `figures_p84v2/stats_20k_seed42.tsv`"),
    "S3_hyperparameters": ("Architectures, hyperparameters and parameter counts",
                           "Layers, widths, heads, dropout, optimizer, learning rate, weight decay, epochs, early stopping, seeds and parameter counts of TCR-IMAP (both variants), Struct-TCS, Seq-TCS, TEIM-Res and AngleNet",
                           "`scripts/19`, `21`, `24b`, `28b`, `53_benchmark_p84v2`, `20`, `36b`, `16b`; parameter counts by instantiation"),
    "S4_protocol_folds": ("Partitioning protocols: fold composition",
                          "Per protocol and fold: n training, n test, clusters on each side, maximum and mean train/test similarity, duplicate families split; allele/intersection coverage limitation",
                          "`protocols.json`, `complex_metadata.tsv`, `tcr_identity_matrix.npy`, `leakage_audit_dual.tsv`, `percomplex_maxsim.tsv`, `validation_fix/B2_split_dryrun.log`"),
    "S5_identity_sweep": ("Identity-threshold sweep: composition and results",
                          "Per threshold: n test, n training per fold, achieved maximum and mean similarity, the `mean_max_sim` upper bound, Δ, CI, p",
                          "`similarity_curve_p84v2.csv`, `similarity_curve_composition{,_folds}_p84v2.csv`"),
    "S6_md_selection_qc": ("Molecular dynamics: system selection and per-replicate quality control",
                           "The four systems (selection basis, allele, peptide, resolution, p84v2 rank) and, per replicate, steps, frames, LINCS warnings, temperature, pressure, density, backbone RMSD, TCR–pMHC contact retention",
                           "`md/MD_SYSTEM_SELECTION_CRITERIA.md`, `percomplex_tcr70_paired_p84v2.tsv`, `characterisation_fix/S2_qc_all12.log`, `md_characterisation_*.csv`, `md_occupancy_docked_frames.csv`"),
    "S7_md_occupancy": ("Molecular dynamics: occupancy statistics",
                        "Pooled, per-system, per-stratum and per-replicate Spearman ρ for TCR-IMAP and the prior with counts; top-k enrichment; crystal-positive occupancy; docked-frame analysis; score–score correlation",
                        "`md_occupancy_{rho,strata,topk,check1_positive,docked_rho}.csv`, `md_score_score_spearman.csv`"),
    "S8_dataset_descriptors": ("Dataset descriptors",
                               "Exact twins, twin pairs differing in peptide, deduplication (132 complexes, median NN 0.967), helix share 88.3 %, β-sheet floor, each with its computation",
                               "`tcr_identity_matrix.npy`, `complex_metadata.tsv`, `identity_bound/handoff6_statements.tsv`, `gdomain182_region_loop_summary_p84v2.csv`"),
}
rows = []
for stem, (title, cap, src) in caps.items():
    csvs = sorted(p.name for p in OUT.glob(f"{stem.split('_')[0]}_*.csv"))
    rows.append([f"[{stem.split('_')[0]}]({stem}.md)", title, cap, src, ", ".join(f"`{c}`" for c in csvs)])
pend = "\n".join(f"- {p}" for p in PENDING) or "- none"
index = f"""# Supporting Information tables: index

Built by `build_si_tables.py` (this directory) from committed tables, reports and stored per-complex arrays. No
training, inference, scoring or GPU. Recomputed statistics use the existing pooled clu70 cluster bootstrap
(`src/viz/figstats.pooled_cluster_bootstrap`, 20,000 resamples, seed 42) on stored per-complex arrays and are asserted
equal to the committed values; the build stops on any mismatch. Log: `build.log`. Numbering follows the task order
(S1–S8), not yet citation order; these tables are distinct from Figures S1–S9.

{md_table(["Table", "Title", "Caption (short)", "Source files", "CSV"], rows)}

## Pending

{pend}
"""
(OUT / "INDEX.md").write_text(index)
print("wrote INDEX.md")
print("PENDING:"); print(pend)
