#!/usr/bin/env python3
"""p84v2 rebuild of the 180-position descriptive audit (40d): TCR x MHC G-domain contacts and the
DSSP region classification, with residues keyed by (IMGT number, insertion code).

40d's definitions are kept: same 174 complexes in 40d's order, same structures, heavy atoms,
4.5 A, same DSSP method (HEADER-prepended temp copy, mkdssp, helix = H/G/I, sheet = E/B) and the
same thresholds (MIN_RESOLVED_FOR_CALL, MIN_MAJORITY_FRAC, imported from 40d itself).
Validation first: the same code with number-only keying (last residue in file wins; DSSP
out[resnum] = ss in key order) must reproduce 40d's stored contact_freq_gdomain180_descriptive.npz
(Y, S, F) and gdomain180_region_classification.csv exactly. Only then are p84v2 outputs written.
p84v2: TCR rows = the 75 p84v2 rows; MHC columns = IMGT 1-90 and 1001-1090 with insertion residues
as their own columns (1061A, 1072A) -> 182.
Outputs (refuse to overwrite): results/tables/contact_freq_gdomain182_descriptive_p84v2.npz,
gdomain182_region_classification_p84v2.csv, gdomain182_region_loop_summary_p84v2.csv
usage: python scripts/53_gdomain182_descriptive_p84v2.py
"""
import importlib.util
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, DSSP
from scipy.spatial.distance import cdist

sys.path.insert(0, ".")
from src.data import p84v2 as P  # noqa: E402

spec = importlib.util.spec_from_file_location("d40", "scripts/40d_gdomain180_descriptive.py")
d40 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d40)      # constants only; main() not called
HELIX, SHEET, MIN_RES, MIN_MAJ, CUT = d40.HELIX_CODES, d40.SHEET_CODES, d40.MIN_RESOLVED_FOR_CALL, d40.MIN_MAJORITY_FRAC, d40.CUTOFF
POS180 = d40.POS180
OUT_NPZ = Path("results/tables/contact_freq_gdomain182_descriptive_p84v2.npz")
OUT_REG = Path("results/tables/gdomain182_region_classification_p84v2.csv")
OUT_SUM = Path("results/tables/gdomain182_region_loop_summary_p84v2.csv")
for f in (OUT_NPZ, OUT_REG, OUT_SUM):
    assert not f.exists(), f"refusing to overwrite existing {f}"

ref = np.load("results/tables/contact_freq_gdomain180_descriptive.npz", allow_pickle=True)
ref_reg = pd.read_csv("results/tables/gdomain180_region_classification.csv")
PDBS = [str(p) for p in ref["pdb"]]                                   # 40d's own order
TL82 = [str(t) for t in ref["tcr_position_labels"]]
lab = P.load(); rows84 = list(zip([s.split(":")[0] for s in lab["row_label"]], lab["row_imgt_number"], lab["row_insertion_code"]))
R84 = {k: i for i, k in enumerate(rows84)}
idx = pd.read_csv("data/interim/stcrdab_index.tsv", sep="\t").drop_duplicates("pdb").set_index("pdb")
parser = PDBParser(QUIET=True)


def dssp(pdb, chain):
    src = Path(f"data/structures/stcrdab/imgt/{pdb}.pdb").read_text()
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as tf:
        tf.write("HEADER    SYNTHETIC HEADER FOR DSSP PARSING\n"); tf.write(src); tmp = tf.name
    try:
        model = parser.get_structure(pdb, tmp)[0]
        dd = DSSP(model, tmp, dssp="mkdssp")
        keys = [k for k in dd.keys() if k[0] == chain]
        return [(k[1][1], k[1][2].strip(), dd[k][2]) for k in keys]      # in DSSP key order
    except Exception:
        return None
    finally:
        Path(tmp).unlink(missing_ok=True)


def classify(keys, codes):
    out = []
    for key in keys:
        num = key if isinstance(key, int) else key[0]
        cs = codes.get(key, []); n = len(cs); cc = Counter("helix" if c in HELIX else "sheet" if c in SHEET else "other" for c in cs)
        if n < MIN_RES:
            region, frac = "unresolved/uncertain", np.nan
        else:
            top, tn = cc.most_common(1)[0]; frac = tn / n
            region = ("unresolved/uncertain" if frac <= MIN_MAJ else
                      ("alpha1 helix" if num <= 90 else "alpha2 helix") if top == "helix" else
                      "beta-sheet floor" if top == "sheet" else "other G-domain")
        out.append(dict(imgt_position=key if isinstance(key, int) else f"{key[0]}{key[1]}", domain="alpha1" if num <= 90 else "alpha2",
                        n_resolved=n, n_helix=cc["helix"], n_sheet=cc["sheet"], n_other=cc["other"],
                        majority_frac=round(float(frac), 3) if not np.isnan(frac) else "", region=region))
    return pd.DataFrame(out)


def contacts(tres_list, mres_list):
    """residue lists -> boolean contact matrix (min heavy-atom distance <= CUT)."""
    ta, to, ma, mo = [], [], [], []
    for i, x in enumerate(tres_list):
        for at in x:
            if at.element != "H": ta.append(at.coord); to.append(i)
    for j, x in enumerate(mres_list):
        for at in x:
            if at.element != "H": ma.append(at.coord); mo.append(j)
    M = np.zeros((len(tres_list), len(mres_list)), bool)
    if ta and ma:
        h = cdist(np.array(ta, float), np.array(ma, float)) <= CUT
        ii, jj = np.nonzero(h); M[np.array(to)[ii], np.array(mo)[jj]] = True
    return M


# ---------------- pass 1: structures, both keyings ----------------
codes82, codes84, colkeys84 = {}, {}, set()
Y82 = np.zeros((174, 84, 180), np.int8); S82 = np.zeros_like(Y82)
store = {}
for n, pdb in enumerate(PDBS):
    r = idx.loc[pdb]; mch = str(r.mhc_chain1)
    model = parser.get_structure(pdb, f"data/structures/stcrdab/imgt/{pdb}.pdb")[0]
    ss = dssp(pdb, mch)
    assert ss is not None, f"DSSP failed for {pdb}"
    last = {}
    for num, ic, code in ss:
        last[num] = code                                                  # 40d: out[resnum] = ss, last wins
        if num in set(POS180):
            codes84.setdefault((num, ic), []).append(code)
    for p in POS180:
        if p in last: codes82.setdefault(p, []).append(last[p])
    mstd = [x for x in model[mch] if x.id[0] == " " and x.id[1] in set(POS180)]
    for x in mstd: colkeys84.add((x.id[1], x.id[2].strip()))
    tstd = {}
    for tag, ch in (("a", r.Achain), ("b", r.Bchain)):
        tstd[tag] = [x for x in model[str(ch)] if x.id[0] == " "]
    store[pdb] = (mstd, tstd)
    # p82v1-mode contacts (40d logic)
    mres = {x.id[1]: x for x in mstd}
    rows_res, rows_idx = [], []
    for tag in ("a", "b"):
        tres = {x.id[1]: x for x in tstd[tag]}
        for seg, (lo, _) in d40.SEG_RANGE.items():
            key = f"{seg}{tag}"; off = sum(d40.SEG_LENS[k] for k in d40.SEG_ORDER[:d40.SEG_ORDER.index(key)])
            for p in range(lo, lo + d40.SEG_LENS[key]):
                if p in tres: rows_res.append(tres[p]); rows_idx.append(off + p - lo)
    cols_idx = [j for j, p in enumerate(POS180) if p in mres]
    Cm = contacts(rows_res, [mres[POS180[j]] for j in cols_idx])
    S82[n][np.ix_(rows_idx, cols_idx)] = 1; Y82[n][np.ix_(rows_idx, cols_idx)] = Cm

# ---------------- validation against 40d ----------------
reg82 = classify(POS180, codes82)
import io as _io   # compare the written representation: 40d's CSV stores '' for unresolved majority_frac, read back as NaN
reg82 = pd.read_csv(_io.StringIO(reg82.to_csv(index=False)))
v_y, v_s = bool((Y82 == ref["Y"]).all()), bool((S82 == ref["S"]).all())
cmp = reg82.merge(ref_reg, on="imgt_position", suffixes=("", "_ref"))
v_reg = all((cmp[c].astype(str) == cmp[f"{c}_ref"].astype(str)).all() for c in ("n_resolved", "n_helix", "n_sheet", "n_other", "majority_frac", "region"))
print(f"validation vs 40d: Y {v_y}, S {v_s}, region table {v_reg} (rows {len(cmp)}/180)")
if not v_reg:
    for c in ("n_resolved", "n_helix", "n_sheet", "n_other", "majority_frac", "region"):
        bad = cmp[cmp[c].astype(str) != cmp[f"{c}_ref"].astype(str)]
        if len(bad):
            print(f"   column {c}: {len(bad)} rows differ, e.g.", bad[["imgt_position", c, f"{c}_ref"]].head(6).values.tolist())
if not (v_y and v_s and v_reg and len(cmp) == 180):
    print("STOP: the reimplementation does not reproduce 40d"); sys.exit(1)

# ---------------- p84v2 ----------------
cols84 = sorted(set((p, "") for p in POS180) | colkeys84, key=lambda k: P.imgt_rank(*k))
C84 = {k: j for j, k in enumerate(cols84)}
Y84 = np.zeros((174, len(rows84), len(cols84)), np.int8); S84 = np.zeros_like(Y84)
for n, pdb in enumerate(PDBS):
    mstd, tstd = store[pdb]
    rr, rres = [], []
    for tag in ("a", "b"):
        for x in tstd[tag]:
            for seg, (lo, _) in d40.SEG_RANGE.items():
                key = f"{seg}{tag}"
                if lo <= x.id[1] < lo + d40.SEG_LENS[key]:
                    rr.append(R84[(key, x.id[1], x.id[2].strip())]); rres.append(x)
    cc = [C84[(x.id[1], x.id[2].strip())] for x in mstd]
    Cm = contacts(rres, mstd)
    S84[n][np.ix_(rr, cc)] = 1; Y84[n][np.ix_(rr, cc)] = Cm
reg84 = classify(cols84, codes84)
F84 = (Y84 * S84).sum(0) / np.maximum(S84.sum(0), 1)
col_labels = np.array([f"{p}{i}" for p, i in cols84])
np.savez_compressed(OUT_NPZ, F=F84, Y=Y84, S=S84, pdb=np.array(PDBS), mhc_positions=col_labels,
                    mhc_imgt_number=np.array([p for p, _ in cols84]), mhc_insertion_code=np.array([i for _, i in cols84]),
                    tcr_position_labels=np.array(lab["row_label"]), cutoff=np.array(CUT),
                    candidate_set=np.array("gdomain182_descriptive_p84v2"), excludes_alpha3=np.array(True),
                    provenance=np.array("scripts/53_gdomain182_descriptive_p84v2.py; 40d definitions, (IMGT number, insertion code) keying; "
                                        "validated: number-only mode reproduces 40d exactly"))
reg84.to_csv(OUT_REG, index=False)
region_of = dict(zip(reg84.imgt_position.astype(str), reg84.region))
loops = {}
for i, l in enumerate(lab["row_label"]):
    loops.setdefault(l.split(":")[0], []).append(i)
rows = []
tot = int(Y84.sum())
for loop, ri in loops.items():
    for region in ["alpha1 helix", "alpha2 helix", "beta-sheet floor", "other G-domain", "unresolved/uncertain"]:
        cj = [j for j, l in enumerate(col_labels) if region_of[l] == region]
        if not cj: continue
        nc = int(Y84[:, ri][:, :, cj].sum()); no = int(S84[:, ri][:, :, cj].sum())
        rows.append(dict(loop=loop, region=region, n_mhc_positions=len(cj), n_contacts=nc, n_pairs_observed=no,
                         frac_of_all_contacts=round(nc / tot, 4), contact_rate=round(nc / no, 4) if no else np.nan))
pd.DataFrame(rows).to_csv(OUT_SUM, index=False)
print(f"wrote {OUT_NPZ} {Y84.shape}, {OUT_REG}, {OUT_SUM}")
print("region counts p84v2:", reg84.region.value_counts().to_dict())
print("region counts p82v1:", ref_reg.region.value_counts().to_dict())
