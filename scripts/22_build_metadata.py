#!/usr/bin/env python3
"""
Phase 1: complex-level metadata, duplicate detection, TCR sequence clustering.

Produces data/processed/complex_metadata.tsv with one row per complex:
  pdb, date, study_id, allele, peptide, TRAV/TRBV, CDR sequences,
  cluster IDs at 70/80/90% identity, duplicate-family ID, resolution.
"""
import sys, re
import numpy as np, pandas as pd
from Levenshtein import ratio
sys.path.insert(0,'.')
from src.models.grip import SEG_ORDER

CDRS = ["cdr1a","cdr2a","cdr3a","cdr1b","cdr2b","cdr3b"]

def single_linkage(seqs, thr):
    """cluster indices by pairwise identity >= thr (single linkage)"""
    n = len(seqs)
    parent = list(range(n))
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    def union(a,b):
        ra,rb=find(a),find(b)
        if ra!=rb: parent[ra]=rb
    for i in range(n):
        for j in range(i+1,n):
            if ratio(seqs[i], seqs[j]) >= thr:
                union(i,j)
    roots = {}
    out = []
    for i in range(n):
        r = find(i)
        if r not in roots: roots[r]=len(roots)
        out.append(roots[r])
    return np.array(out)

def main():
    meta = pd.read_parquet("data/processed/pair_meta.parquet")
    idx = pd.read_csv("data/interim/stcrdab_index.tsv", sep="\t").drop_duplicates("pdb").set_index("pdb")
    rows = []
    for _, m in meta.iterrows():
        r = idx.loc[m.pdb]
        rows.append(dict(
            pdb=m.pdb, allele=m.allele, trav=m.trav, trbv=m.trbv,
            peptide=str(r.antigen_name)[:60],
            date=str(r.date), resolution=r.resolution,
            study_id=str(r.authors)[:80],           # proxy for source publication
            docking_angle=r.docking_angle,
            **{c: m[c] for c in CDRS}))
    df = pd.DataFrame(rows)

    df["tcr_seq"] = df[CDRS].agg("".join, axis=1)
    df["cdr3_pair"] = df.cdr3a + "|" + df.cdr3b

    # duplicate family: same TCR + same allele + same peptide
    df["dup_family"] = df.groupby(["cdr3_pair","allele","peptide"]).ngroup()
    dup_sizes = df.dup_family.value_counts()
    print(f"complexes: {len(df)} | duplicate families: {df.dup_family.nunique()}")
    print(f"  families with >1 structure: {(dup_sizes>1).sum()} covering {dup_sizes[dup_sizes>1].sum()} complexes")

    for thr, name in ((0.70,"clu70"), (0.80,"clu80"), (0.90,"clu90")):
        df[name] = single_linkage(df.tcr_seq.tolist(), thr)
        print(f"  TCR clusters at {int(thr*100)}% identity: {df[name].nunique()}")

    # date parsing for temporal split
    def parse_date(s):
        try: return pd.to_datetime(s, format="%m/%d/%y")
        except Exception: return pd.NaT
    df["dep_date"] = df.date.map(parse_date)
    df.loc[df.dep_date > "2026-12-31", "dep_date"] -= pd.DateOffset(years=100)
    print(f"\ndate range: {df.dep_date.min()} .. {df.dep_date.max()}")
    for cut in ("2018-01-01","2020-01-01","2022-01-01"):
        n = (df.dep_date >= cut).sum()
        print(f"  structures deposited >= {cut}: {n}")

    print(f"\nstudies (author-string proxy): {df.study_id.nunique()}")
    print(f"alleles: {df.allele.nunique()} | TRBV: {df.trbv.nunique()} | TRAV: {df.trav.nunique()}")
    df.to_csv("data/processed/complex_metadata.tsv", sep="\t", index=False)
    print("\nwrote data/processed/complex_metadata.tsv")

if __name__=="__main__":
    main()
