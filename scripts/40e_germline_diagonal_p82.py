"""Germline diagonal under the 180-position DSSP-classified descriptive set.

Descriptive geometry only -- no model, no split, no candidate-set selection on outcome.
Asks one question per TCR loop: does it contact the alpha1 helix or the alpha2 helix more,
aggregated to clu70 clusters so that near-duplicate complexes count once.

Statistic: per-cluster contact-rate difference d_c = rate_alpha1(c) - rate_alpha2(c).
Reported as the MEAN over clusters with a cluster bootstrap CI on that same mean
(2000 resamples, p floored at 1/2000 per project convention), plus a sign test as a
distribution-free secondary. Holm correction across the six loops.

This concerns docking-footprint geometry (which helix a loop sits over). It is NOT the
withdrawn stereotypy claim and does NOT constitute two-step-docking confirmation.
"""
import numpy as np, pandas as pd
from math import comb

N_BOOT, SEED = 2000, 42
rng = np.random.default_rng(SEED)

z = np.load("results/tables/contact_freq_gdomain180_descriptive.npz", allow_pickle=True)
Y, S, pos, tcr, pdb = z["Y"], z["S"], z["mhc_positions"], z["tcr_position_labels"], z["pdb"]
reg = (pd.read_csv("results/tables/gdomain180_region_classification.csv")
         .set_index("imgt_position")["region"].reindex(pos).values)
A1, A2 = reg == "alpha1 helix", reg == "alpha2 helix"
loop_of = np.array([str(t).split(":")[0] for t in tcr])

meta = pd.read_csv("data/processed/complex_metadata.tsv", sep="\t")
kcol = [c for c in meta.columns if c.lower() in ("pdb", "pdb_id", "id")][0]
clu = meta.set_index(kcol)["clu70"].reindex([str(p) for p in pdb]).values
assert not pd.isna(clu).any()
uclu = np.unique(clu)

def per_cluster(li, M):
    k = np.array([Y[clu == c][:, li][:, :, M].sum() for c in uclu], float)
    n = np.array([(S[clu == c][:, li][:, :, M] == 1).sum() for c in uclu], float)
    return k, n

rows = []
for L in ["cdr1a", "cdr2a", "cdr3a", "cdr1b", "cdr2b", "cdr3b"]:
    li = loop_of == L
    k1, n1 = per_cluster(li, A1); k2, n2 = per_cluster(li, A2)
    ok = (n1 > 0) & (n2 > 0)
    d = k1[ok] / n1[ok] - k2[ok] / n2[ok]
    boot = np.array([rng.choice(d, d.size, replace=True).mean() for _ in range(N_BOOT)])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    p_boot = max(2 * min((boot <= 0).mean(), (boot >= 0).mean()), 1 / N_BOOT)
    npos, nneg = int((d > 0).sum()), int((d < 0).sum())
    m = npos + nneg
    p_sign = min(1.0, 2 * sum(comb(m, i) for i in range(max(npos, nneg), m + 1)) / 2**m) if m else 1.0
    rows.append(dict(loop=L,
                     contacts_a1=int(k1.sum()), resolved_a1=int(n1.sum()),
                     contacts_a2=int(k2.sum()), resolved_a2=int(n2.sum()),
                     rate_a1=k1.sum()/n1.sum(), rate_a2=k2.sum()/n2.sum(),
                     clusters_informative=int(ok.sum()), clusters_tied=int((d == 0).sum()),
                     n_a1_gt_a2=npos, n_a2_gt_a1=nneg,
                     mean_diff=d.mean(), ci_lo=lo, ci_hi=hi,
                     p_bootstrap=p_boot, p_sign=p_sign,
                     prefers="alpha1" if d.mean() > 0 else "alpha2"))

df = pd.DataFrame(rows)
# Holm across the six loops, on the bootstrap p
order = df["p_bootstrap"].values.argsort()
holm = np.empty(len(df)); running = 0.0
for rank, i in enumerate(order):
    running = max(running, (len(df) - rank) * df["p_bootstrap"].values[i])
    holm[i] = min(running, 1.0)
df["p_holm"] = holm

df.to_csv("results/tables/germline_diagonal_p82.csv", index=False)
pd.set_option("display.width", 200)
print(df[["loop", "prefers", "contacts_a1", "contacts_a2", "rate_a1", "rate_a2",
          "clusters_informative", "n_a1_gt_a2", "n_a2_gt_a1", "clusters_tied",
          "mean_diff", "ci_lo", "ci_hi", "p_bootstrap", "p_sign", "p_holm"]].to_string(index=False))
print("\nwrote results/tables/germline_diagonal_p82.csv")
print(f"seed={SEED} n_boot={N_BOOT} bootstrap p floored at {1/N_BOOT} -> report as p < 0.001")
