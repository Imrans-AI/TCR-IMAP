"""Per-complex metrics and cluster-bootstrap confidence intervals."""
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, matthews_corrcoef

def per_complex(y, p, mask):
    """y,p,mask: (N, ...) arrays. Returns dict of per-complex metric arrays."""
    out = {k: [] for k in ("auc","auprc","p_at_L","p_at_10","mcc","brier","n_pos")}
    for i in range(len(y)):
        m = mask[i].astype(bool).ravel()
        yi, pi = y[i].ravel()[m], p[i].ravel()[m]
        npos = int(yi.sum())
        out["n_pos"].append(npos)
        if npos == 0 or npos == len(yi):
            for k in ("auc","auprc","p_at_L","p_at_10","mcc","brier"): out[k].append(np.nan)
            continue
        order = np.argsort(-pi)
        out["auc"].append(roc_auc_score(yi, pi))
        out["auprc"].append(average_precision_score(yi, pi))
        out["p_at_L"].append(yi[order[:npos]].mean())
        out["p_at_10"].append(yi[order[:min(10,len(yi))]].mean())
        out["mcc"].append(matthews_corrcoef(yi, (pi >= 0.5).astype(int)))
        out["brier"].append(np.mean((pi - yi)**2))
    return {k: np.array(v, dtype=float) for k, v in out.items()}

def cluster_bootstrap_diff(m_a, m_b, clusters, n_boot=2000, seed=0):
    """Paired cluster bootstrap of mean(m_a) - mean(m_b). Returns (diff, lo, hi, p)."""
    rng = np.random.default_rng(seed)
    ok = ~(np.isnan(m_a) | np.isnan(m_b))
    a, b, c = m_a[ok], m_b[ok], np.asarray(clusters)[ok]
    uc = np.unique(c)
    obs = a.mean() - b.mean()
    diffs = np.empty(n_boot)
    for i in range(n_boot):
        pick = rng.choice(uc, len(uc), replace=True)
        idx = np.concatenate([np.where(c == u)[0] for u in pick])
        diffs[i] = a[idx].mean() - b[idx].mean()
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    p = 2 * min((diffs <= 0).mean(), (diffs >= 0).mean())
    return obs, lo, hi, max(p, 1/n_boot)

def paired_permutation(m_a, m_b, clusters, n_perm=5000, seed=0):
    rng = np.random.default_rng(seed)
    ok = ~(np.isnan(m_a) | np.isnan(m_b))
    d = (m_a - m_b)[ok]; c = np.asarray(clusters)[ok]
    obs = d.mean()
    uc = np.unique(c)
    cnt = 0
    for _ in range(n_perm):
        flip = {u: rng.choice([-1,1]) for u in uc}
        s = np.array([flip[x] for x in c])
        if abs((d*s).mean()) >= abs(obs): cnt += 1
    return obs, (cnt+1)/(n_perm+1)
