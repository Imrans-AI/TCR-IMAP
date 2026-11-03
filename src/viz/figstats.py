"""Figure-time statistics shared by the p82v1 figure scripts (Figs 2, 4, 6, SI peptide).

Nothing here trains, predicts or scores. Every function operates on per-complex arrays that
already exist as artifacts.

pooled_cluster_bootstrap() is the estimator of src.eval.metrics.cluster_bootstrap_diff
(resample clu70 clusters with replacement, pool their complexes, mean paired difference,
percentile CI, two-sided p floored at 1/n_boot), vectorised so 20,000 resamples are cheap.
For the same (n_boot, seed) it consumes the RNG identically and returns identical values;
the figure scripts assert this against the original at 2,000 resamples before using it.
"""
import re
import zlib
import numpy as np

N_BOOT = 20000
SEED = 42


def _clean(a, b, clusters):
    a = np.asarray(a, float); b = np.asarray(b, float); c = np.asarray(clusters)
    ok = ~(np.isnan(a) | np.isnan(b))
    return a[ok], b[ok], c[ok]


def _cluster_counts(n_clusters, n_boot, seed):
    """(n_boot, K) multiplicity of each cluster per resample; same draws as calling
    rng.choice(uc, K, replace=True) once per resample."""
    rng = np.random.default_rng(seed)
    cnt = np.empty((n_boot, n_clusters), np.int32)
    for i in range(n_boot):
        pick = rng.choice(n_clusters, n_clusters, replace=True)
        cnt[i] = np.bincount(pick, minlength=n_clusters)
    return cnt


def _p_from(diffs, n_boot):
    p = 2 * min((diffs <= 0).mean(), (diffs >= 0).mean())
    return max(p, 1 / n_boot)


def pooled_cluster_bootstrap(a, b, clusters, n_boot=N_BOOT, seed=SEED, return_draws=False):
    """mean(a) - mean(b) over complexes, clusters resampled and pooled.
    Returns dict(diff, lo, hi, p, n, k, n_boot, seed); with return_draws=True also 'draws',
    the n_boot resampled differences, so further quantiles (e.g. simultaneous intervals) come
    from the SAME resamples as the 95 % CI and p."""
    a, b, c = _clean(a, b, clusters)
    uc, inv = np.unique(c, return_inverse=True)
    k = len(uc)
    sa = np.bincount(inv, weights=a, minlength=k)
    sb = np.bincount(inv, weights=b, minlength=k)
    nn = np.bincount(inv, minlength=k).astype(float)
    cnt = _cluster_counts(k, n_boot, seed).astype(float)
    tot = cnt @ nn
    diffs = (cnt @ sa) / tot - (cnt @ sb) / tot
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    out = dict(diff=float(a.mean() - b.mean()), lo=float(lo), hi=float(hi),
               p=float(_p_from(diffs, n_boot)), n=int(len(a)), k=int(k),
               n_boot=n_boot, seed=seed)
    if return_draws:
        out["draws"] = diffs
    return out


def cluster_macro_bootstrap(a, b, clusters, n_boot=N_BOOT, seed=SEED):
    """Cluster-macro sensitivity: per-cluster mean of each method (clusters equally
    weighted), difference of the cluster-level means; clusters resampled with replacement.
    Returns dict(diff, macro_a, macro_b, lo, hi, p, n, k, n_boot, seed)."""
    a, b, c = _clean(a, b, clusters)
    uc, inv = np.unique(c, return_inverse=True)
    k = len(uc)
    nn = np.bincount(inv, minlength=k).astype(float)
    ma = np.bincount(inv, weights=a, minlength=k) / nn
    mb = np.bincount(inv, weights=b, minlength=k) / nn
    d = ma - mb
    cnt = _cluster_counts(k, n_boot, seed).astype(float)
    diffs = (cnt @ d) / k
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return dict(diff=float(d.mean()), macro_a=float(ma.mean()), macro_b=float(mb.mean()),
                lo=float(lo), hi=float(hi), p=float(_p_from(diffs, n_boot)),
                n=int(len(a)), k=int(k), n_boot=n_boot, seed=seed)


def holm(pvals):
    """Holm step-down adjusted p-values, same order as input, capped at 1."""
    p = np.asarray(pvals, float)
    m = len(p)
    order = np.argsort(p)
    adj = np.empty(m)
    run = 0.0
    for r, i in enumerate(order):
        run = max(run, (m - r) * p[i])
        adj[i] = min(run, 1.0)
    return adj


def ci_excludes_zero(lo, hi):
    return not (lo <= 0 <= hi)


def parse_cell(s):
    """'+0.0390 [+0.0158,+0.0626] p=0.0010' -> (diff, lo, hi, p)."""
    m = re.match(r"\s*([+-]?\d*\.\d+)\s*\[([+-]?\d*\.\d+),\s*([+-]?\d*\.\d+)\]\s*p=([\d.]+)", s)
    if not m:
        raise ValueError(f"unparseable comparison cell: {s!r}")
    return tuple(float(x) for x in m.groups())


def pdf_font_check(path):
    """Return (n_type3, has_arialmt, basefonts) by scanning the PDF's object streams."""
    raw = open(path, "rb").read()
    chunks = [raw]
    for m in re.finditer(rb"stream\r?\n(.*?)\r?\nendstream", raw, re.S):
        try:
            chunks.append(zlib.decompress(m.group(1)))
        except Exception:
            pass
    blob = b"\n".join(chunks)
    n_t3 = len(re.findall(rb"/Subtype\s*/Type3", blob))
    bases = sorted({x.decode("latin1") for x in re.findall(rb"/BaseFont\s*/([A-Za-z0-9+\-_,]+)", blob)})
    return n_t3, any("ArialMT" in b for b in bases), bases
