#!/usr/bin/env python3
"""Spearman correlation between the two scores of Figure 7 - TCR-IMAP (pred_model) and the pair-frequency prior
(pred_freq) - over the crystal-negative resolved pairs of the four MD systems, per system and pooled.

Recorded because the Figure 7 legend cites it ("the two scores are only moderately correlated, Spearman 0.33-0.40")
and that value had been computed ad hoc and never written down. Scores only: no occupancy is read. The score file is
checked against the MD run manifest's scores_sha256, so these are exactly the scores the run used; the pair universe
is the analysis's own (stats.neg_mask on the p84v2 labels, also checked by sha256). Descriptive: no p value (pairs
within a system are not independent, and nothing is tested here).

writes results/tables/md_score_score_spearman.csv
usage: python scripts/56_score_score_md.py
"""
import hashlib
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, "md/analysis")
import stats as st                           # the analysis's own crystal-negative mask

T = "results/tables/"
SYSTEMS = ["5yxn", "7n1f", "3vxs", "8eo8"]
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
man = pd.read_csv(T + "md_occupancy_manifest.csv").set_index("item")["value"]
SCORES, LABELS = T + "pairmaps_tcr70_p84v2.npz", "data/processed/pair_maps_p84v2.npz"
assert sha(SCORES) == man["scores_sha256"] and sha(LABELS) == man["labels_sha256"]
pm, lab = np.load(SCORES, allow_pickle=True), np.load(LABELS, allow_pickle=True)
lpdb = [str(x) for x in lab["pdb"]]
assert [str(x) for x in pm["pdb"]] == lpdb

rows, pooled = [], ([], [])
rho_tab = pd.read_csv(T + "md_occupancy_rho.csv", float_precision="round_trip")
for p in SYSTEMS:
    i = lpdb.index(p)
    Y, S = lab["Y"][i], lab["S"][i]
    assert (pm["Y"][i] == Y).all() and (pm["S"][i] == S).all()
    m = st.neg_mask(Y, S)
    a, b = pm["pred_model"][i][m].astype(float), pm["pred_freq"][i][m].astype(float)
    n_tab = int(rho_tab[(rho_tab.level == "system") & (rho_tab.pdb == p) & (rho_tab.score == "TCR-IMAP")].n_crystal_negative.iloc[0])
    assert len(a) == n_tab, (p, len(a), n_tab)              # the same pair universe as the occupancy tables
    rows.append(dict(level="system", pdb=p, n_crystal_negative=len(a), n_distinct_tcrimap=len(np.unique(a)),
                     n_distinct_prior=len(np.unique(b)), spearman_rho=float(spearmanr(a, b).correlation)))
    pooled[0].append(a); pooled[1].append(b)
a, b = np.concatenate(pooled[0]), np.concatenate(pooled[1])
rows.append(dict(level="pooled", pdb="all four", n_crystal_negative=len(a), n_distinct_tcrimap=len(np.unique(a)),
                 n_distinct_prior=len(np.unique(b)), spearman_rho=float(spearmanr(a, b).correlation)))
D = pd.DataFrame(rows)
D.to_csv(T + "md_score_score_spearman.csv", index=False)
print(D.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
sysr = D[D.level == "system"].spearman_rho
print(f"\nper-system range {sysr.min():.2f}-{sysr.max():.2f}; pooled {D.spearman_rho.iloc[-1]:.2f}  "
      f"(legend had 'Spearman 0.33-0.40', source unrecorded)")
