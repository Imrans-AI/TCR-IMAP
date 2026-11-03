#!/usr/bin/env python3
"""Part A of the v3 revision: verify the external review's findings against the frozen
artifacts and the runners' own definitions. Read-only; writes one TSV per finding to
results/review/v3_partA/ and prints the tables that appear in V3_PART_A_VERIFICATION.md.

Training membership is taken as the runners define it -- every non-test complex
(scripts/28b_benchmark_p82v1.py line 77, `rest = np.setdiff1d(arange(n), te)`) -- not as the
protocol builder's group assignment suggests. That distinction is the whole of finding A1.

usage: python scripts/59_verify_review_findings.py
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from Levenshtein import ratio, distance

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results/review/v3_partA"; OUT.mkdir(parents=True, exist_ok=True)
MIN_FOLD = 3          # 28b's `if len(te)<3: continue` guard

df = pd.read_csv(ROOT/"data/processed/complex_metadata.tsv", sep="\t")
proto = json.load(open(ROOT/"data/splits_v2/protocols.json"))
npz = np.load(ROOT/"data/processed/pair_maps_p84v2.npz", allow_pickle=True)
assert list(npz["pdb"]) == list(df.pdb.values), "pair_maps_p84v2 is not in metadata order"
df["mhc_window"] = ["".join(r) for r in npz["mhc_aa"]]


def a1_a2():
    """A1: the double protocol's held-out groups stay in training. A2: the allele protocol
    holds out labels, not MHC windows."""
    rows = []
    for pname in ("double", "allele"):
        for fi, fold in enumerate(proto[pname]):
            te = df[df.pdb.isin(fold)]; tr = df[~df.pdb.isin(fold)]
            if len(te) == 0: continue
            tr_alleles = set(tr.allele.dropna())
            tr_windows = set(tr.mhc_window)
            tr_windows_unannotated = set(tr[tr.allele.isna()].mhc_window)
            shares_cluster = te.clu70.isin(set(tr.clu70))
            shares_allele = te.allele.notna() & te.allele.isin(tr_alleles)
            rows.append(dict(
                protocol=pname, fold=fi, evaluated=len(te) >= MIN_FOLD,
                n_test=len(te), n_train=len(tr), test_clusters=te.clu70.nunique(),
                shares_training_cluster=int(shares_cluster.sum()),
                shares_known_training_allele=int(shares_allele.sum()),
                unseen_in_both=int(((~shares_cluster) & (~shares_allele)).sum()),
                exact_window_in_training=int(te.mhc_window.isin(tr_windows).sum()),
                exact_window_from_unannotated=int(te.mhc_window.isin(tr_windows_unannotated).sum())))
    t = pd.DataFrame(rows)
    t.to_csv(OUT/"a1_a2_protocol_overlap.tsv", sep="\t", index=False)
    print("\n=== A1/A2: overlap between test folds and the training set they are scored against ===")
    print(t.to_string(index=False))
    ev = t[t.evaluated].groupby("protocol")[
        ["n_test", "shares_training_cluster", "shares_known_training_allele",
         "unseen_in_both", "exact_window_in_training", "exact_window_from_unannotated"]].sum()
    print("\ntotals over evaluated folds (n >= 3):"); print(ev.to_string())
    print(f"\ndistinct 84-residue MHC windows over the 174 complexes: {df.mhc_window.nunique()}")
    return t


def cluster_disjoint_validation(clu, rest, rng, frac=0.15):
    """Verbatim from scripts/53_benchmark_p84v2.py VALIDATION_HELPER."""
    cl = clu[rest]; uc = np.unique(cl); target = frac*len(rest); sel, n = [], 0
    for c in uc[rng.permutation(len(uc))]:
        if n >= target: break
        sel.append(c); n += int((cl == c).sum())
    return rest[np.isin(cl, sel)]


def a3(protocols=("tcr70",)):
    """A3: cluster disjointness of the primary split against training AND validation.
    Reproduces the runner's draw: one default_rng(0) per protocol, one draw per fold in fold
    order, folds below MIN_FOLD skipped before the draw."""
    clu, seqs = df.clu70.values, df.tcr_seq.values
    pdb2i = {p: i for i, p in enumerate(df.pdb)}
    rows = []
    for pname in protocols:
        rng = np.random.default_rng(0)
        for fi, fold in enumerate(proto[pname]):
            te = np.array([pdb2i[p] for p in fold if p in pdb2i])
            if len(te) < MIN_FOLD: continue
            rest = np.setdiff1d(np.arange(len(df)), te)
            va = cluster_disjoint_validation(clu, rest, rng, 0.15)
            tr = np.setdiff1d(rest, va)
            rows.append(dict(
                protocol=pname, fold=fi, n_test=len(te), n_nontest=len(rest),
                n_train=len(tr), n_val=len(va),
                val_pct_of_nontest=round(100*len(va)/len(rest), 1),
                val_clusters=len(set(clu[va])),
                overlap_test_train=len(set(clu[tr]) & set(clu[te])),
                overlap_test_val=len(set(clu[va]) & set(clu[te])),
                overlap_test_nontest=len(set(clu[rest]) & set(clu[te])),
                overlap_train_val=len(set(clu[tr]) & set(clu[va])),
                max_sim_test_to_train=round(max(max(ratio(seqs[i], seqs[k]) for k in tr) for i in te), 4),
                max_sim_test_to_nontest=round(max(max(ratio(seqs[i], seqs[k]) for k in rest) for i in te), 4)))
    t = pd.DataFrame(rows)
    t.to_csv(OUT/"a3_tcr70_disjointness.tsv", sep="\t", index=False)
    print("\n=== A3: primary split, against actual training and validation membership ===")
    print(t.to_string(index=False))
    return t


def a4():
    """A4: what Levenshtein.ratio actually normalises. The equation form is right; the
    distance is misnamed. ratio() is the Indel normalised similarity (substitution costs 2)."""
    def indel(a, b):
        la, lb = len(a), len(b)
        dp = [[0]*(lb+1) for _ in range(la+1)]
        for i in range(1, la+1):
            for j in range(1, lb+1):
                dp[i][j] = dp[i-1][j-1]+1 if a[i-1] == b[j-1] else max(dp[i-1][j], dp[i][j-1])
        return la + lb - 2*dp[la][lb]
    rows = []
    for a, b in [("ABC","ABD"), ("CAT","COT"), ("ABCDE","ABXDE"), ("KITTEN","SITTING"), ("AAAA","AAAB")]:
        n = len(a) + len(b)
        rows.append(dict(a=a, b=b, ratio=round(ratio(a, b), 4),
                         levenshtein_form=round((n - distance(a, b))/n, 4),
                         indel_form=round((n - indel(a, b))/n, 4)))
    t = pd.DataFrame(rows)
    t["matches_indel"] = np.isclose(t.ratio, t.indel_form)
    t["matches_levenshtein"] = np.isclose(t.ratio, t.levenshtein_form)
    t.to_csv(OUT/"a4_ratio_semantics.tsv", sep="\t", index=False)
    print("\n=== A4: which distance does Levenshtein.ratio normalise? ===")
    print(t.to_string(index=False))
    print(f"matches Indel in {int(t.matches_indel.sum())}/{len(t)} cases, "
          f"Levenshtein in {int(t.matches_levenshtein.sum())}/{len(t)}")
    return t


def a5():
    """A5: receptor-species composition. The cohort is human by MHC, not by receptor:
    scripts/02_build_index.py filters mhc_chain1_organism only."""
    ix = pd.read_csv(ROOT/"data/interim/stcrdab_index.tsv", sep="\t")
    m = ix[ix.pdb.isin(df.pdb)].drop_duplicates("pdb").set_index("pdb")
    assert len(m) == len(df), f"matched {len(m)} of {len(df)}"
    norm = lambda s: m[s].astype(str).str.strip().str.lower()
    print("\n=== A5: species composition of the 174 complexes ===")
    for c in ("alpha_organism", "beta_organism", "mhc_chain1_organism", "mhc_chain2_organism", "mhc_type"):
        print(f"\n{c}:"); print(norm(c).value_counts().to_string())
    nonhuman = m[(norm("alpha_organism") != "homo sapiens") | (norm("beta_organism") != "homo sapiens")]
    print("\ncomplexes with a non-human TCR chain:")
    print(nonhuman[["alpha_organism", "beta_organism", "mhc_chain1_organism", "mhc_type"]].to_string())
    nonhuman.to_csv(OUT/"a5_nonhuman_receptors.tsv", sep="\t")
    return nonhuman


if __name__ == "__main__":
    a1_a2(); a3(); a4(); a5()
    print(f"\nwrote {OUT}")
