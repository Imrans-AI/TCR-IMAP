#!/usr/bin/env python3
"""Compare two executions of the p84v2 queue, table by table. No training, no inference, no re-scoring.

OLD = the backup directory a `--force` run of scripts/53_run_all_p84v2.sh made of the previous execution:
      <old>/<stage>/<output file>, <old>/<stage>/<cache dir>/, <old>/derived_snapshot/<derived table>.
NEW = a directory laid out like results/tables (default results/tables): <new>/<output file>,
      <new>/_partial/p84v2/<cache dir>/.

The table list is read from the queue script's STAGES (the 20 stage outputs: 9 CSV, 1 JSONL, 12 .npz), plus the
three tables derived outside the queue (pairwise_ours_p84v2.npz, the 13th git-ignored .npz;
percomplex_tcr70_paired_p84v2.tsv; cluster_macro_p84v2.csv), plus every file in every stage cache.

For every file: byte-identical or not. Where not, every value is compared:
  * CSV/TSV cells: numbers directly; text cells by their numeric tokens when the surrounding text is the same
    (e.g. "+0.0504 [+0.0193,+0.0823] p=0.0020"), otherwise as a text difference;
  * JSONL: record by record, key by key;
  * .npz: array by array, element by element (NaN in the same place counts as equal).
Each differing value is classed as NOISE (|difference| <= 1e-6) or REAL (> 1e-6); text, structure (shape, keys,
columns) and NaN-placement differences are REAL. Every differing value is written to diffs.tsv.gz with its
magnitude (for an array with more than --max-rows differing elements, the largest --max-rows are written and the
rest are summarised in keys.tsv; --all writes every one).

Headline check, at the precision the paper reports (macro-AUPRC, differences and CI bounds to 4 decimals, p as
style.pstr() prints it, docking-angle MAE to 2 and r / R^2 to 3 decimals, win counts, showcase complexes, Holm
survivors). Two kinds: the RECORDED values the tables hold (2,000 resamples), and the paper-standard 20,000-
resample (seed 42) statistics recomputed here from the stored per-complex arrays of each execution (per-residue
benchmark, pairwise vs prior and TEIM-Res, cluster-macro, ablations with Holm), exactly as
scripts/55_stats_p84v2.py does. The pairwise rows need the derived tables re-derived from the NEW run first
(scripts/55_stats_p84v2.py); if the NEW derived tables are older than the NEW queue outputs they are reported STALE.

Outputs in --out: summary.md, files.tsv, keys.tsv, diffs.tsv.gz, headline.tsv.
Exit status: 0 every value equal or within 1e-6 and no headline change; 1 real differences but no headline
change; 2 a headline value changes at paper precision; 3 a file is missing or stale.

  python scripts/55_compare_p84v2_runs.py --old results/review/train_p84v2_logs/forced_backup_<time> \
      --new results/tables --out results/review/clean_rerun/compare_<time>
"""
import argparse, gzip, hashlib, io, json, math, os, re, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from viz.figstats import pooled_cluster_bootstrap as pcb, cluster_macro_bootstrap as cmb, holm  # noqa: E402
from viz.style import pstr  # noqa: E402

NOISE = 1e-6
QUEUE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "53_run_all_p84v2.sh")
DERIVED = ["pairwise_ours_p84v2.npz", "percomplex_tcr70_paired_p84v2.tsv", "cluster_macro_p84v2.csv"]
NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")
PROTOS = ["random", "double", "allele", "temporal", "tcr70"]
METHODS = ["Struct-TCS", "Seq-TCS", "PosFreq", "NN"]
ABL = ["no_mhc", "no_tcr", "shuffle_tcr", "shuffle_mhc", "no_germline", "no_cdr3", "no_cross", "dot_head"]


def stages_from_queue(path=QUEUE):
    """[(stage, [output basenames], [cache dir names])] parsed from the queue's STAGES block."""
    s = open(path).read()
    block = s[s.index("STAGES=(\n"):]
    block = block[:block.index("\n)\n")]
    out = []
    for line in re.findall(r'^"([^"]+)"$', block, flags=re.M):
        name, _runner, outputs, caches = line.split("|")
        out.append((name, [os.path.basename(o) for o in outputs.split()], [os.path.basename(c) for c in caches.split()]))
    assert len(out) == 8, out
    return out


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class Diffs:
    """Collects every differing value; per-file and per-key summaries."""
    def __init__(self, max_rows, write_all):
        self.rows, self.keys, self.max_rows, self.all = [], [], max_rows, write_all

    def numeric(self, f, key, idx, old, new, label=None):
        """old/new: aligned 1-D float arrays (NaN allowed); idx: labels (list) or a function i -> label."""
        lab = idx if callable(idx) else (lambda i: idx[i])
        old = np.asarray(old, float); new = np.asarray(new, float)
        both_nan = np.isnan(old) & np.isnan(new)
        nan_mis = np.isnan(old) ^ np.isnan(new)
        d = np.where(both_nan | nan_mis, 0.0, np.abs(new - old))
        diff = (d > 0) | nan_mis
        n_diff, n_real = int(diff.sum()), int(((d > NOISE) | nan_mis).sum())
        rec = dict(file=f, key=key, n=int(old.size), n_diff=n_diff, n_real=n_real, nan_mismatch=int(nan_mis.sum()),
                   max_abs=float(d.max()) if d.size else 0.0, mean_abs_of_diffs=float(d[diff].mean()) if n_diff else 0.0,
                   where_max="" if not n_diff else str(lab(int(np.argmax(np.where(diff, d, -1))))),
                   cls="equal" if n_diff == 0 else ("REAL" if n_real else "NOISE"), written=0)
        if n_diff:
            sel = np.where(diff)[0]
            if not self.all and len(sel) > self.max_rows:
                sel = sel[np.argsort(-d[sel], kind="stable")[: self.max_rows]]
            for i in sel:
                self.rows.append((f, key, str(lab(int(i))), repr(float(old[i])), repr(float(new[i])),
                                  "nan-mismatch" if nan_mis[i] else f"{d[i]:.3e}",
                                  "REAL" if (nan_mis[i] or d[i] > NOISE) else "NOISE"))
            rec["written"] = int(len(sel))
        self.keys.append(rec)
        return rec

    def text(self, f, key, where, old, new):
        self.rows.append((f, key, str(where), repr(old)[:300], repr(new)[:300], "text", "REAL"))

    def structure(self, f, key, msg):
        self.keys.append(dict(file=f, key=key, n=0, n_diff=1, n_real=1, nan_mismatch=0, max_abs=float("nan"),
                              mean_abs_of_diffs=float("nan"), where_max=msg, cls="STRUCTURE", written=1))
        self.rows.append((f, key, "-", "-", "-", msg, "REAL"))


def tokens(s):
    return NUM.sub("#", s), [float(x) for x in NUM.findall(s)]


def as_float(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def compare_table(f, po, pn, D):
    sep = "\t" if f.endswith((".tsv", ".tab")) else ","
    a = pd.read_csv(po, sep=sep, dtype=str, keep_default_na=False)
    b = pd.read_csv(pn, sep=sep, dtype=str, keep_default_na=False)
    if list(a.columns) != list(b.columns) or a.shape != b.shape:
        D.structure(f, "*", f"columns/shape differ: {list(a.columns)} {a.shape} vs {list(b.columns)} {b.shape}")
        return
    rk = a.columns[0]
    texts = 0
    for c in a.columns:
        idx, ov, nv = [], [], []
        for i in range(len(a)):
            x, y = a[c].iat[i], b[c].iat[i]
            where = f"{rk}={a[rk].iat[i]}"
            fx, fy = as_float(x), as_float(y)
            if fx is not None and fy is not None:
                idx.append(where); ov.append(fx); nv.append(fy)
            elif x == y:
                continue
            else:
                (sx, tx), (sy, ty) = tokens(x), tokens(y)
                if sx == sy and len(tx) == len(ty):
                    for t, (u, v) in enumerate(zip(tx, ty)):
                        idx.append(f"{where} token{t}"); ov.append(u); nv.append(v)
                else:
                    D.text(f, c, where, x, y); texts += 1
        if idx:
            D.numeric(f, c, idx, ov, nv)
    if texts:
        D.keys.append(dict(file=f, key="(text cells)", n=0, n_diff=texts, n_real=texts, nan_mismatch=0,
                           max_abs=float("nan"), mean_abs_of_diffs=float("nan"), where_max="see diffs", cls="TEXT",
                           written=texts))


def flat(obj, pre=""):
    if isinstance(obj, dict):
        for k in sorted(obj):
            yield from flat(obj[k], f"{pre}.{k}" if pre else str(k))
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            yield from flat(v, f"{pre}[{i}]")
    else:
        yield pre, obj


def compare_jsonl(f, po, pn, D):
    A = [json.loads(l) for l in open(po) if l.strip()]
    B = [json.loads(l) for l in open(pn) if l.strip()]
    if len(A) != len(B):
        D.structure(f, "*", f"{len(A)} vs {len(B)} records"); return
    idx, ov, nv = [], [], []
    for r, (x, y) in enumerate(zip(A, B)):
        fx, fy = dict(flat(x)), dict(flat(y))
        if set(fx) != set(fy):
            D.structure(f, f"record {r}", f"keys differ: {sorted(set(fx) ^ set(fy))}"); continue
        lab = "/".join(str(x.get(k)) for k in ("split_by", "seed") if k in x) or str(r)
        for k in fx:
            u, v = fx[k], fy[k]
            if isinstance(u, (int, float)) and isinstance(v, (int, float)) and not isinstance(u, bool):
                idx.append(f"{lab}:{k}"); ov.append(u); nv.append(v)
            elif u != v:
                D.text(f, k, lab, u, v)
    D.numeric(f, "(all numeric fields)", idx, ov, nv)


def compare_npz(f, po, pn, D):
    A, B = np.load(po, allow_pickle=True), np.load(pn, allow_pickle=True)
    if sorted(A.files) != sorted(B.files):
        D.structure(f, "*", f"keys differ: {sorted(set(A.files) ^ set(B.files))}")
    for k in sorted(set(A.files) & set(B.files)):
        x, y = A[k], B[k]
        if x.shape != y.shape:
            D.structure(f, k, f"shape {x.shape} vs {y.shape}"); continue
        if x.dtype.kind in "fiub" and y.dtype.kind in "fiub":
            shp = x.shape
            D.numeric(f, k, (lambda i, shp=shp: str(tuple(int(v) for v in np.unravel_index(i, shp))) if shp else "()"),
                      x.reshape(-1).astype(float), y.reshape(-1).astype(float))
        else:
            xs, ys = x.reshape(-1), y.reshape(-1)
            neq = [i for i in range(xs.size) if xs[i] != ys[i]]
            for i in neq[:1000]:
                D.text(f, k, np.unravel_index(i, x.shape) if x.ndim else "()", xs[i], ys[i])
            D.keys.append(dict(file=f, key=k, n=int(xs.size), n_diff=len(neq), n_real=len(neq), nan_mismatch=0,
                               max_abs=float("nan"), mean_abs_of_diffs=float("nan"),
                               where_max="" if not neq else "text", cls="equal" if not neq else "TEXT",
                               written=min(len(neq), 1000)))


def compare_file(f, po, pn, D):
    if not os.path.exists(po) or not os.path.exists(pn):
        return dict(file=f, status="MISSING", old=os.path.exists(po), new=os.path.exists(pn))
    so, sn = sha(po), sha(pn)
    st = dict(file=f, old_sha256=so[:16], new_sha256=sn[:16], size_old=os.path.getsize(po), size_new=os.path.getsize(pn))
    if so == sn:
        st["status"] = "byte-identical"; return st
    k0, r0 = len(D.keys), len(D.rows)
    if f.endswith(".npz"):
        compare_npz(f, po, pn, D)
    elif f.endswith(".jsonl"):
        compare_jsonl(f, po, pn, D)
    else:
        compare_table(f, po, pn, D)
    ks = D.keys[k0:]
    real = any(k["cls"] in ("REAL", "TEXT", "STRUCTURE") for k in ks)
    noise = any(k["cls"] == "NOISE" for k in ks)
    st["status"] = "REAL differences" if real else ("noise only (<= 1e-6)" if noise else "values equal (bytes differ)")
    st["max_abs"] = max([k["max_abs"] for k in ks if isinstance(k["max_abs"], float) and not math.isnan(k["max_abs"])] or [0.0])
    st["n_diff_values"] = sum(k["n_diff"] for k in ks)
    st["classes"] = ",".join(sorted({k["cls"] for k in ks if k["cls"] != "equal"})) or "equal"
    return st


# ------------------------------------------------------------------ headline values at paper precision
def f4(x): return f"{x:+.4f}" if x is not None else "NA"


def cell(s):
    """'+0.0504 [+0.0193,+0.0823] p=0.0020' -> paper-precision text."""
    v = [float(t) for t in NUM.findall(s)]
    return f"{v[0]:+.4f} [{v[1]:+.4f}, {v[2]:+.4f}] {pstr(v[3])}"


def headline(side):
    """side: dict basename -> path (stage outputs, caches dir map 'cache:<name>' -> dir, derived)."""
    H = {}
    g = lambda b: side.get(b)
    if g("benchmark_p84v2.csv"):
        b = pd.read_csv(g("benchmark_p84v2.csv")).set_index("protocol")
        for p in b.index:
            for m in ["Ours"] + METHODS:
                H[f"recorded per-residue {p} {m} macro-AUPRC"] = f"{b.loc[p, m]:.4f}"
            for m in METHODS:
                H[f"recorded per-residue {p} TCR-IMAP vs {m}"] = cell(b.loc[p, f"vs_{m}"])
    if g("protocol_main_p84v2.csv"):
        b = pd.read_csv(g("protocol_main_p84v2.csv")).set_index("protocol")
        for p in b.index:
            r = b.loc[p]
            H[f"recorded pairwise {p} TCR-IMAP / prior"] = f"{r.model_auprc:.4f} / {r.freq_auprc:.4f}"
            H[f"recorded pairwise {p} vs prior"] = f"{r['diff']:+.4f} [{r.ci_lo:+.4f}, {r.ci_hi:+.4f}] {pstr(r.pval)}"
    if g("ablations_p84v2.csv"):
        b = pd.read_csv(g("ablations_p84v2.csv")).set_index("ablation")
        for a in b.index:
            r = b.loc[a]
            H[f"recorded ablation {a}"] = f"{r.auprc:.4f}; {r['diff']:+.4f} [{r.lo:+.4f}, {r.hi:+.4f}] {pstr(r.p)}"
    if g("peptide_aware_p84v2.csv"):
        b = pd.read_csv(g("peptide_aware_p84v2.csv")).set_index("scope")
        for s in b.index:
            r = b.loc[s]
            H[f"recorded peptide {s}"] = f"{r.no_pep:.4f} -> {r.pep:.4f}; {r['diff']:+.4f} [{r.lo:+.4f}, {r.hi:+.4f}] {pstr(r.pval)}"
    if g("similarity_curve_p84v2.csv"):
        b = pd.read_csv(g("similarity_curve_p84v2.csv"))
        for _, r in b.iterrows():
            H[f"recorded identity sweep < {r.threshold:g}"] = (f"{r.model:.4f} vs {r.freq:.4f}; {r['diff']:+.4f} "
                                                               f"[{r.lo:+.4f}, {r.hi:+.4f}] {pstr(r.p)}")
    for tag in ("4p5A", "6A"):
        if g(f"structcs_setting_p84v2_{tag}.csv"):
            b = pd.read_csv(g(f"structcs_setting_p84v2_{tag}.csv")).set_index("setting")
            for s in b.index:
                r = b.loc[s]
                H[f"recorded comparator subset {tag} {s}"] = (f"n={r.n_test}; {r['Ours_auprc']:.4f} / {r['Struct-TCS_auprc']:.4f}; "
                                                              f"{cell(r['vs_Struct-TCS'])}")
    if g("angle_runs_p84v2.jsonl"):
        A = pd.DataFrame([json.loads(l) for l in open(g("angle_runs_p84v2.jsonl")) if l.strip()])
        for s, d in A.groupby("split_by"):
            H[f"docking angle {s} (MAE / r / R2 ranges, n)"] = (f"{d.mae.min():.2f}-{d.mae.max():.2f} / {d.pearson.min():.3f}-"
                                                              f"{d.pearson.max():.3f} / {d.r2.min():+.3f} to {d.r2.max():+.3f}, n={d.n.iloc[0]}")
    if g("teim_res_pairwise_p84v2.npz"):
        H["TEIM-Res tcr70 pairwise macro-AUPRC"] = f"{np.nanmean(np.load(g('teim_res_pairwise_p84v2.npz'), allow_pickle=True)['auprc']):.4f}"
    if g("cluster_macro_p84v2.csv"):
        b = pd.read_csv(g("cluster_macro_p84v2.csv"))
        for _, r in b.iterrows():
            if pd.notna(r.p_cluster_macro):
                H[f"recorded cluster-macro {r.task} {r.model}"] = (f"{r.paired_diff_vs_GRIP_cluster_macro:+.4f} [{r.ci_lo_cluster_macro:+.4f}, "
                                                                 f"{r.ci_hi_cluster_macro:+.4f}] {pstr(r.p_cluster_macro)}")

    # ---- paper-standard 20,000 resamples, seed 42, from the stored per-complex arrays
    def s20(r): return f"{r['diff']:+.4f} [{r['lo']:+.4f}, {r['hi']:+.4f}] {pstr(r['p'])}"
    for p in PROTOS:
        if g(f"percomplex_{p}_p84v2.npz"):
            z = np.load(g(f"percomplex_{p}_p84v2.npz"), allow_pickle=True)
            for m in METHODS:
                H[f"20k per-residue {p} TCR-IMAP vs {m}"] = s20(pcb(z["Ours_auprc"], z[f"{m}_auprc"], z["clu70"]))
            if p == "tcr70":
                for m in ("Struct-TCS", "PosFreq"):
                    H[f"20k cluster-macro per-residue tcr70 vs {m}"] = s20(cmb(z["Ours_auprc"], z[f"{m}_auprc"], z["clu70"]))
    if g("pairwise_ours_p84v2.npz") and g("percomplex_tcr70_paired_p84v2.tsv"):
        po = np.load(g("pairwise_ours_p84v2.npz"), allow_pickle=True)
        pr = pd.read_csv(g("percomplex_tcr70_paired_p84v2.tsv"), sep="\t")
        H["20k pairwise tcr70 vs prior"] = s20(pcb(po["auprc"], pr.prior.values, po["clu70"]))
        H["20k cluster-macro pairwise tcr70 vs prior"] = s20(cmb(po["auprc"], pr.prior.values, po["clu70"]))
        H["pairwise tcr70 TCR-IMAP / prior macro-AUPRC"] = f"{np.nanmean(po['auprc']):.4f} / {pr.prior.mean():.4f}"
        H["pairwise tcr70 wins over prior"] = f"{int((pr.delta > 0).sum())}/{len(pr)}"
        o = pr.sort_values("delta", ascending=False, kind="stable")
        H["Fig 6 showcase (max / min delta)"] = f"{o.pdb.iloc[0]} {o.delta.iloc[0]:+.3f} / {o.pdb.iloc[-1]} {o.delta.iloc[-1]:+.3f}"
        if g("teim_res_pairwise_p84v2.npz"):
            te = np.load(g("teim_res_pairwise_p84v2.npz"), allow_pickle=True)
            H["20k pairwise tcr70 vs TEIM-Res"] = s20(pcb(po["auprc"], te["auprc"], po["clu70"]))
            H["pairwise tcr70 wins over TEIM-Res"] = f"{int((po['auprc'] > te['auprc']).sum())}/{len(po['auprc'])}"
    cache = g("cache:29b_ablations")
    if cache and os.path.exists(os.path.join(cache, "full.npz")) and g("pairwise_ours_p84v2.npz"):
        full = np.load(os.path.join(cache, "full.npz"))
        # cluster labels in the 29b cache order (idx 0..173 = pairwise_ours order, as 41b asserts)
        assert (full["idx"] == np.arange(174)).all()
        cl = np.load(g("pairwise_ours_p84v2.npz"), allow_pickle=True)["clu70"] if g("pairwise_ours_p84v2.npz") else None
        res = {a: pcb(full["auprc"], np.load(os.path.join(cache, f"{a}.npz"))["auprc"], cl) for a in ABL}
        ph = holm([res[a]["p"] for a in ABL])
        for a, h in zip(ABL, ph):
            H[f"20k ablation {a} (Holm 8)"] = f"{s20(res[a])}; Holm {pstr(h)}"
        H["ablations surviving Holm (20k)"] = ", ".join(a for a, h in zip(ABL, ph) if h < 0.05)
    return H


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--old", required=True, help="forced_backup_<time> directory of the previous execution")
    ap.add_argument("--new", default="results/tables", help="directory laid out like results/tables")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-rows", type=int, default=100000)
    ap.add_argument("--all", action="store_true", help="write every differing element")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    D = Diffs(a.max_rows, a.all)
    files, oldmap, newmap = [], {}, {}
    for stage, outs, caches in stages_from_queue():
        for o in outs:
            oldmap[o], newmap[o] = os.path.join(a.old, stage, o), os.path.join(a.new, o)
            files.append((stage, o, oldmap[o], newmap[o]))
        for c in caches:
            od, nd = os.path.join(a.old, stage, c), os.path.join(a.new, "_partial", "p84v2", c)
            oldmap[f"cache:{c}"], newmap[f"cache:{c}"] = od, nd
            names = sorted(set(os.listdir(od) if os.path.isdir(od) else []) | set(os.listdir(nd) if os.path.isdir(nd) else []))
            if not names:
                files.append((stage, f"{c}/ (cache)", od, nd))
            for n in names:
                files.append((stage, f"_partial/p84v2/{c}/{n}", os.path.join(od, n), os.path.join(nd, n)))
    for d in DERIVED:
        oldmap[d], newmap[d] = os.path.join(a.old, "derived_snapshot", d), os.path.join(a.new, d)
        files.append(("derived", d, oldmap[d], newmap[d]))

    # derived tables must be re-derived from the NEW outputs before they mean anything
    new_outs = [newmap[o] for _, outs, _ in stages_from_queue() for o in outs if os.path.exists(newmap[o])]
    t_new = min(os.path.getmtime(p) for p in new_outs) if new_outs else None
    stale = [d for d in DERIVED if os.path.exists(newmap[d]) and t_new and os.path.getmtime(newmap[d]) < t_new]

    status = []
    for stage, f, po, pn in files:
        st = compare_file(f, po, pn, D)
        st["stage"] = stage
        if stage == "derived" and f in stale:
            st["status"] = "STALE (not re-derived from the new run) / " + st["status"]
        status.append(st)
    S = pd.DataFrame(status)[lambda d: ["stage", "file", "status"] + [c for c in d.columns if c not in ("stage", "file", "status")]]
    S.to_csv(os.path.join(a.out, "files.tsv"), sep="\t", index=False)
    K = pd.DataFrame(D.keys)
    K.to_csv(os.path.join(a.out, "keys.tsv"), sep="\t", index=False)
    with gzip.open(os.path.join(a.out, "diffs.tsv.gz"), "wt") as fh:
        fh.write("file\tkey\tlocation\told\tnew\tabs_diff\tclass\n")
        for r in D.rows:
            fh.write("\t".join(r) + "\n")

    HO, HN = headline(oldmap), headline(newmap)
    hl = [dict(item=k, old=HO.get(k, "absent"), new=HN.get(k, "absent"), changed=HO.get(k) != HN.get(k))
          for k in list(HO) + [k for k in HN if k not in HO]]
    HL = pd.DataFrame(hl)
    HL.to_csv(os.path.join(a.out, "headline.tsv"), sep="\t", index=False)

    n_real = int((K.cls.isin(["REAL", "TEXT", "STRUCTURE"])).sum()) if len(K) else 0
    missing = S[S.status == "MISSING"]
    changed = HL[HL.changed]
    lines = [f"# p84v2 run comparison", "", f"OLD: `{a.old}`  ", f"NEW: `{a.new}`", "",
             "## Files", "", S.status.value_counts().rename_axis("status").reset_index().to_string(index=False), ""]
    for _, r in S[S.status != "byte-identical"].iterrows():
        num = pd.notna(r.get("n_diff_values")) and ("REAL" in str(r.get("classes")) or "NOISE" in str(r.get("classes")))
        lines.append(f"- `{r.file}` ({r.stage}): {r.status}" +
                     (f"; classes {r.classes}; {int(r.n_diff_values)} differing values" if pd.notna(r.get("n_diff_values")) else "") +
                     (f"; max numeric |diff| {r.max_abs:.3e}" if num else ""))
    if stale:
        lines += ["", f"**STALE derived tables** (older than the new queue outputs; re-derive with scripts/55_stats_p84v2.py "
                      f"and compare again): {', '.join(stale)}"]
    lines += ["", "## Differences by array / column", ""]
    nz = K[K.n_diff > 0] if len(K) else K
    lines += ["none" if not len(nz) else nz[["file", "key", "n", "n_diff", "n_real", "nan_mismatch", "max_abs", "where_max", "cls", "written"]].to_string(index=False)]
    lines += ["", f"Every differing value is in diffs.tsv.gz ({len(D.rows)} rows; NOISE = |diff| <= {NOISE:g}, REAL otherwise)."]
    lines += ["", "## Headline values at the precision the paper reports", ""]
    if changed.empty:
        lines.append(f"**No headline value changes at paper precision** ({len(HL)} values compared: recorded 2,000-resample "
                     f"tables, 20,000-resample statistics from the per-complex arrays, docking-angle ranges, win counts, "
                     f"showcase complexes, Holm survivors).")
    else:
        lines.append(f"**{len(changed)} of {len(HL)} headline values CHANGE at paper precision:**")
        lines += [f"- {r['item']}: `{r['old']}` -> `{r['new']}`" for _, r in changed.iterrows()]
    if stale:
        lines.append("(pairwise, win-count, showcase and ablation-cluster rows use the STALE derived tables listed above)")
    verdict = 3 if (len(missing) or stale) else 2 if len(changed) else 1 if n_real else 0
    lines += ["", f"Exit status {verdict}: " + {0: "all values equal or within 1e-6; no headline change",
                                               1: "real differences (> 1e-6) but no headline change at paper precision",
                                               2: "headline values change at paper precision",
                                               3: "missing or stale files"}[verdict]]
    open(os.path.join(a.out, "summary.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return verdict


if __name__ == "__main__":
    sys.exit(main())
