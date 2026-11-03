#!/usr/bin/env python3
"""Descriptive-only 84x180 pooled contact-frequency map over the full MHC G-domain
(alpha1 IMGT 1-90 + alpha2 IMGT 1001-1090), for the 174-complex primary dataset (6uz1
excluded, matching every other p82v1 artifact).

This is NOT a predictive analysis: no model is trained, no train/test split is used, no
comparator is invoked. It computes one pooled, descriptive quantity -- observed contact
frequency per (TCR position, MHC position) cell across all 174 complexes -- purely to
test whether real contact density exists outside the alpha1/alpha2 helix faces, i.e. the
question the main audit flagged as "not supported as a discovery claim" on the biased
28-position candidate set.

SCOPE, STATED EXPLICITLY: the 180-position G-domain set covers only the alpha1+alpha2
platform domains. It EXCLUDES the alpha3 immunoglobulin-like domain (IMGT >= 1091) --
the same domain where the 6uz1 diagnostic found the majority of that (excluded) complex's
real contacts. Nothing in this script's output can support any claim about the alpha3
domain, in either direction.

Region classification of the 180 candidate positions uses a structure-derived, documented
method (not a memorized numbering convention): DSSP secondary structure, computed directly
from each of the 174 crystal structures' MHC chain, majority-voted per IMGT position across
all structures in which that position is resolved. See classify_positions() docstring for
the exact rule and thresholds.
"""
import subprocess, sys, tempfile, time
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np, pandas as pd
from Bio.PDB import PDBParser, DSSP
from Bio.SeqUtils import seq1
sys.path.insert(0, ".")
from src.models.grip import SEG_ORDER, SEG_LENS, TCR_LEN

SEG_RANGE = {"cdr1": (27, 38), "cdr2": (56, 65), "cdr3": (105, 117)}
POS180 = list(range(1, 91)) + list(range(1001, 1091))
EXCLUDE_PDBS = {"6uz1"}
CUTOFF = 4.5
MIN_RESOLVED_FOR_CALL = 30      # documented threshold: need >=30/174 structures resolving a
                                  # position before classifying its secondary structure at all
MIN_MAJORITY_FRAC = 0.50         # documented threshold: winning SS-category must exceed 50%
                                  # of resolved calls, else "unresolved/uncertain"

HELIX_CODES = set("HGI")
SHEET_CODES = set("EB")
OTHER_CODES = set("TSP-  ")      # turn, bend, polyproline, coil/none


def dssp_for_structure(pdb, mhc_chain):
    """Run DSSP on one structure's MHC chain via a HEADER-prepended temp copy (the raw
    STCRDab IMGT-renumbered files have no HEADER record, which both mkdssp and Biopython's
    format sniffer require; verified necessary and sufficient in this session before
    writing this function)."""
    raw = Path(f"data/structures/stcrdab/imgt/{pdb}.pdb")
    with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as tf:
        tf.write("HEADER    SYNTHETIC HEADER FOR DSSP PARSING\n")
        tf.write(raw.read_text())
        tmp_path = tf.name
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb, tmp_path)
        model = structure[0]
        dssp = DSSP(model, tmp_path, dssp="mkdssp")
    except Exception as e:
        Path(tmp_path).unlink(missing_ok=True)
        return None
    out = {}
    for key in dssp.keys():
        if key[0] != mhc_chain:
            continue
        resnum = key[1][1]
        ss = dssp[key][2]
        out[resnum] = ss
    Path(tmp_path).unlink(missing_ok=True)
    return out


def classify_positions(idx_df):
    """For every position in POS180, collect the DSSP secondary-structure code observed at
    that IMGT position across every one of the 174 primary-set structures where it is
    resolved, then assign one of five categories:

      alpha1 helix       -- majority code in {H,G,I}, position in 1-90, majority frac > 0.5,
                             resolved in >= MIN_RESOLVED_FOR_CALL structures
      alpha2 helix       -- same, position in 1001-1090
      beta-sheet floor   -- majority code in {E,B} (either domain), same thresholds
      other G-domain     -- majority code is turn/bend/coil/other, same thresholds
      unresolved/uncertain -- resolved in < MIN_RESOLVED_FOR_CALL structures, OR no
                             category reaches > MIN_MAJORITY_FRAC of resolved calls
    """
    per_pos_codes = defaultdict(list)
    n_structs = 0
    t0 = time.time()
    for _, r in idx_df.iterrows():
        ss = dssp_for_structure(r.pdb, str(r.mhc_chain1))
        if ss is None:
            print(f"  DSSP failed for {r.pdb}, skipping this structure for classification")
            continue
        n_structs += 1
        for p in POS180:
            if p in ss:
                per_pos_codes[p].append(ss[p])
    print(f"DSSP ran successfully on {n_structs}/{len(idx_df)} structures ({time.time()-t0:.1f}s)")

    rows = []
    for p in POS180:
        codes = per_pos_codes.get(p, [])
        n_resolved = len(codes)
        cat_counts = Counter()
        for c in codes:
            if c in HELIX_CODES: cat_counts["helix"] += 1
            elif c in SHEET_CODES: cat_counts["sheet"] += 1
            else: cat_counts["other"] += 1
        if n_resolved < MIN_RESOLVED_FOR_CALL:
            region = "unresolved/uncertain"
            majority_frac = np.nan
        else:
            top_cat, top_n = cat_counts.most_common(1)[0]
            majority_frac = top_n / n_resolved
            if majority_frac <= MIN_MAJORITY_FRAC:
                region = "unresolved/uncertain"
            elif top_cat == "helix":
                region = "alpha1 helix" if p <= 90 else "alpha2 helix"
            elif top_cat == "sheet":
                region = "beta-sheet floor"
            else:
                region = "other G-domain"
        rows.append(dict(imgt_position=p, domain=("alpha1" if p <= 90 else "alpha2"),
                          n_resolved=n_resolved, n_helix=cat_counts["helix"],
                          n_sheet=cat_counts["sheet"], n_other=cat_counts["other"],
                          majority_frac=round(float(majority_frac), 3) if not np.isnan(majority_frac) else "",
                          region=region))
    return pd.DataFrame(rows), n_structs


def build_contact_map(idx_df):
    """Same per-atom-distance build logic already used throughout this project's pipeline
    (18b_build_pair_labels_p82v1.py, 26_broader_mhc.py), applied to the 180-position
    G-domain set. Purely geometric -- no model, no training, no split."""
    parser = PDBParser(QUIET=True)
    Y, S, pdbs = [], [], []
    for _, r in idx_df.iterrows():
        f = Path(f"data/structures/stcrdab/imgt/{r.pdb}.pdb")
        if not f.exists(): continue
        try: model = parser.get_structure(r.pdb, str(f))[0]
        except Exception: continue
        try: mchain = model[str(r.mhc_chain1)]
        except KeyError: continue
        mres = {x.id[1]: x for x in mchain if x.id[0] == " "}
        if sum(p in mres for p in POS180) < len(POS180) - 8:  # allow a bit more slack at 180 width
            continue
        y = np.zeros((TCR_LEN, len(POS180)), np.int8)
        s = np.zeros((TCR_LEN, len(POS180)), np.int8)
        ok = True
        for tag, ch in (("a", r.Achain), ("b", r.Bchain)):
            try: tchain = model[str(ch)]
            except KeyError: ok = False; break
            tres = {x.id[1]: x for x in tchain if x.id[0] == " "}
            for seg, (lo, _) in SEG_RANGE.items():
                key = f"{seg}{tag}"
                off = sum(SEG_LENS[k] for k in SEG_ORDER[:SEG_ORDER.index(key)])
                for p in range(lo, lo + SEG_LENS[key]):
                    if p not in tres: continue
                    for jj, mp in enumerate(POS180):
                        if mp not in mres: continue
                        s[off + (p - lo), jj] = 1
                        hit = any(np.linalg.norm(at.coord - bt.coord) <= CUTOFF
                                  for at in tres[p] if at.element != "H"
                                  for bt in mres[mp] if bt.element != "H")
                        if hit: y[off + (p - lo), jj] = 1
        if not ok or y.sum() == 0:
            continue
        Y.append(y); S.append(s); pdbs.append(r.pdb)
    return np.stack(Y), np.stack(S), pdbs


def main():
    t_all = time.time()
    out_npz = "results/tables/contact_freq_gdomain180_descriptive.npz"
    out_regions = "results/tables/gdomain180_region_classification.csv"
    out_summary = "results/tables/gdomain180_region_loop_summary.csv"
    for p in (out_npz, out_regions, out_summary):
        assert not Path(p).exists(), f"refusing to overwrite existing {p}"

    idx = pd.read_csv("data/interim/stcrdab_index.tsv", sep="\t")
    primary174 = set(pd.read_csv("data/processed/complex_metadata.tsv", sep="\t").pdb)
    idx = idx[idx.pdb.isin(primary174)].reset_index(drop=True)
    assert len(idx) == 174 and "6uz1" not in idx.pdb.values
    print(f"primary dataset: {len(idx)} complexes (6uz1 excluded, matches p82v1 throughout)")

    print("\n--- structural region classification (DSSP-based) ---")
    region_df, n_structs_dssp = classify_positions(idx)
    region_df.to_csv(out_regions, index=False)
    print(f"wrote {out_regions}")
    print(region_df.region.value_counts().to_string())

    # cross-check against this project's existing, hardcoded "helix" convention
    # (IMGT alpha1 50-90 + alpha2 1050-1090, used throughout 18b/26/27)
    existing_helix_positions = set(range(50, 91)) | set(range(1050, 1091))
    dssp_helix_positions = set(region_df[region_df.region.isin(["alpha1 helix", "alpha2 helix"])].imgt_position)
    print(f"\nexisting hardcoded helix window: {len(existing_helix_positions)} positions")
    print(f"DSSP-classified helix positions (this analysis): {len(dssp_helix_positions)} positions")
    print(f"  in both: {len(existing_helix_positions & dssp_helix_positions)}")
    print(f"  in existing window but NOT DSSP-helix: {sorted(existing_helix_positions - dssp_helix_positions)}")
    print(f"  DSSP-helix but OUTSIDE existing window: {sorted(dssp_helix_positions - existing_helix_positions)}")

    print("\n--- building 84x180 contact map (geometric, no model) ---")
    Y, S, pdbs = build_contact_map(idx)
    print(f"complexes retained for the contact map: {len(pdbs)} / {len(idx)}")
    assert len(pdbs) == 174, f"expected 174 complexes in the contact map, got {len(pdbs)}"
    assert Y.shape == S.shape == (174, 84, 180)

    F = (Y.astype(float) * S.astype(float)).sum(0) / np.maximum(S.astype(float).sum(0), 1)
    pos_rate = (Y.astype(float) * S.astype(float)).sum() / S.astype(float).sum()
    print(f"Y/S shape {Y.shape}, overall positive rate {pos_rate:.6f}")

    tcr_labels = []
    for tag in ("a", "b"):
        for seg, (lo, _) in SEG_RANGE.items():
            key = f"{seg}{tag}"
            for p in range(lo, lo + SEG_LENS[key]):
                tcr_labels.append(f"{key}:{p}")

    np.savez_compressed(out_npz, F=F, Y=Y, S=S, pdb=np.array(pdbs),
                         mhc_positions=np.array(POS180, dtype=np.int64),
                         tcr_position_labels=np.array(tcr_labels),
                         cutoff=np.array(CUTOFF), candidate_set=np.array("gdomain180_descriptive"),
                         n_structs_dssp=np.array(n_structs_dssp),
                         excludes_alpha3=np.array(True), alpha3_definition=np.array("IMGT >= 1091"))
    print(f"wrote {out_npz}")

    # ---- quantify contacts / mass by region and by CDR loop ----
    region_of = dict(zip(region_df.imgt_position, region_df.region))
    regions_ordered = ["alpha1 helix", "alpha2 helix", "beta-sheet floor", "other G-domain", "unresolved/uncertain"]
    loop_bounds = {}
    off = 0
    for tag in ("a", "b"):
        for seg, _ in SEG_RANGE.items():
            key = f"{seg}{tag}"
            loop_bounds[key] = (off, off + SEG_LENS[key]); off += SEG_LENS[key]

    rows = []
    total_contacts = int(Y.sum())
    total_pairs_observed = int(S.sum())
    for loop, (lo, hi) in loop_bounds.items():
        for region in regions_ordered:
            cols = [jj for jj, p in enumerate(POS180) if region_of.get(p) == region]
            if not cols: continue
            y_sub = Y[:, lo:hi, :][:, :, cols]
            s_sub = S[:, lo:hi, :][:, :, cols]
            n_contact = int(y_sub.sum())
            n_observed = int(s_sub.sum())
            rows.append(dict(loop=loop, region=region, n_mhc_positions=len(cols),
                              n_contacts=n_contact, n_pairs_observed=n_observed,
                              frac_of_all_contacts=round(n_contact/total_contacts, 4) if total_contacts else 0.0,
                              contact_rate=round(n_contact/n_observed, 4) if n_observed else np.nan))
    summary = pd.DataFrame(rows)
    summary.to_csv(out_summary, index=False)
    print(f"\nwrote {out_summary}")
    print(summary.to_string(index=False))

    print("\n--- region totals across all loops ---")
    region_totals = summary.groupby("region").agg(
        n_mhc_positions=("n_mhc_positions", "first"),
        n_contacts=("n_contacts", "sum"), n_pairs_observed=("n_pairs_observed", "sum")).reindex(regions_ordered)
    region_totals["frac_of_all_contacts"] = (region_totals.n_contacts / total_contacts).round(4)
    region_totals["contact_rate"] = (region_totals.n_contacts / region_totals.n_pairs_observed).round(4)
    print(region_totals.to_string())
    print(f"\ntotal observed contacts across the whole 84x180 map: {total_contacts}")
    print(f"total observed (resolved) pairs: {total_pairs_observed}")
    print(f"\nEXPLICIT SCOPE STATEMENT: this candidate set (IMGT 1-90, 1001-1090) excludes the "
          f"alpha3 domain (IMGT >= 1091) entirely. No cell above reflects alpha3 occupancy or "
          f"contact, and no conclusion here extends to that region.")
    print(f"\ndone in {time.time()-t_all:.1f}s")

if __name__ == "__main__":
    main()
