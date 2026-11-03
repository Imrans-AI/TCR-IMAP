import re, glob, pandas as pd
rows = [pd.read_csv(f, sep="\t") for f in glob.glob("data/structures/stcrdab/summary/*.tsv")]
df = pd.concat(rows, ignore_index=True)
print("raw rows:", len(df), "| unique PDBs:", df.pdb.nunique())

m = (df.TCRtype.eq("abTCR") & df.antigen_type.eq("peptide") &
     df.mhc_type.eq("MH1") & df.mhc_chain1_organism.eq("homo sapiens"))
f = df[m].drop_duplicates(subset="pdb", keep="first").copy()
print("human MH1 peptide abTCR, deduped:", len(f))

pat = re.compile(r"HLA[- ]?([ABC])\*?\s?(\d{2}):?(\d{2})?", re.I)
def allele(s):
    mm = pat.search(str(s))
    return f"HLA-{mm.group(1).upper()}*{mm.group(2)}:{mm.group(3) or '01'}" if mm else None
f["allele"] = f.compound.map(allele)

print("\nallele parsed:", f.allele.notna().sum(), "/", len(f))
print(f.allele.value_counts().head(12))
print("\nTRBV:", f.beta_subgroup.nunique(), "| TRAV:", f.alpha_subgroup.nunique())
print(f.beta_subgroup.value_counts().head(8))
f.to_csv("data/interim/stcrdab_index.tsv", sep="\t", index=False)
