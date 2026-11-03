"""Residue mapping for the MD occupancy analysis on the p84v2 axes: each of the 75 TCR rows and 84 MHC
columns -> one residue of the CHARMM-GUI topology, by global sequence alignment, never by number.

Difference from mdmap.py (which is the validated p82v1-layout module and is NOT edited): the crystal side
is keyed by **(IMGT number, insertion code)** instead of by number alone, and the axes are the p84v2 ones.
The crystal residue of each row and column is not re-derived here at all - it is read from the label build
itself (`data/processed/pair_maps_p84v2.npz`, arrays `tcr_resid` / `mhc_resid`, which store the residue id
string, e.g. "112A" or "1061A", that the build assigned to that row/column in that complex). The residue's
one-letter code is then checked against the build's stored token (`tcr_aa` / `mhc_aa`), so a map that
disagreed with the labels could not be built silently.

Everything else - per-chain global alignment, splitting MD residues wherever an atom name repeats (CHARMM-GUI
drops insertion codes, so crystal Ala1061 and Ala1061A both become resid 1061 and MDAnalysis merges them),
and the rigid-superposition coordinate verification of every mapped residue - is imported unchanged from
mdmap.py.

Reads no trajectory and no score.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.SeqUtils import seq1

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mdmap import (ROOT, SYSTEMS, COORD_TOL, HIS, align, kabsch, md_residues, chain_segids,  # noqa: E402
                   chain_roles, crystal_model, md_protein_universe)

LABELS = ROOT / "data/processed/pair_maps_p84v2.npz"
_Z = np.load(LABELS, allow_pickle=True)
LPDB = [str(p) for p in _Z["pdb"]]
ROW_LABEL = [str(x) for x in _Z["row_label"]]
ROW_SEG = [str(x) for x in _Z["row_segment"]]
COL_LABEL = [str(x) for x in _Z["col_label"]]
N_ROW, N_COL = len(ROW_LABEL), len(COL_LABEL)
assert (N_ROW, N_COL) == (75, 84), (N_ROW, N_COL)
assert float(_Z["cutoff"]) == 4.5


def slot_labels_p84v2():
    """The 75 TCR row labels ('cdr1a:27', 'cdr3b:112A', ...) and the 84 MHC column labels ('50', '1061A', ...)."""
    return list(ROW_LABEL), list(COL_LABEL)


def crystal_slots_p84v2(pdb, model, roles):
    """(tcr, mhc): the Bio residue of each p84v2 row/column in this complex, or None.

    Taken from the label build's own per-complex residue ids and cross-checked against its tokens."""
    i = LPDB.index(pdb)
    tres_id = [str(x) for x in _Z["tcr_resid"][i]]
    mres_id = [str(x) for x in _Z["mhc_resid"][i]]
    tres_aa = [str(x) for x in _Z["tcr_aa"][i]]
    mres_aa = [str(x) for x in _Z["mhc_aa"][i]]

    def index_chain(ch):
        return {f"{r.id[1]}{r.id[2].strip()}": r for r in model[ch] if r.id[0] == " "}

    chains = {t: index_chain(roles[f"TCR{t}"]) for t in "ab"}
    mchain = index_chain(roles["MHC"])

    tcr = []
    for k, rid in enumerate(tres_id):
        if not rid:
            tcr.append(None); continue
        t = ROW_SEG[k][-1]                      # 'cdr3b' -> 'b'
        r = chains[t].get(rid)
        assert r is not None, f"{pdb} row {ROW_LABEL[k]}: residue {rid} not in chain {roles[f'TCR{t}']}"
        assert seq1(r.get_resname()) == tres_aa[k], (pdb, ROW_LABEL[k], rid, r.get_resname(), tres_aa[k])
        tcr.append(r)
    mhc = []
    for j, rid in enumerate(mres_id):
        if not rid:
            mhc.append(None); continue
        r = mchain.get(rid)
        assert r is not None, f"{pdb} column {COL_LABEL[j]}: residue {rid} not in chain {roles['MHC']}"
        assert seq1(r.get_resname()) == mres_aa[j], (pdb, COL_LABEL[j], rid, r.get_resname(), mres_aa[j])
        mhc.append(r)
    return tcr, mhc


def build_map_p84v2(pdb):
    """Residue map + G1 facts for one system on the p84v2 axes.
    Returns (map DataFrame, per-chain report, chain->segid, role->chain, MD residues by segid, fit rmsd)."""
    model = crystal_model(pdb)
    roles = chain_roles(pdb)
    segids = chain_segids(pdb, model)
    u = md_protein_universe(pdb)
    X = u.atoms.positions
    tcr, mhc = crystal_slots_p84v2(pdb, model, roles)
    tlabs, mlabs = slot_labels_p84v2()
    fmt = lambda r: f"{r.get_resname()}{r.id[1]}{r.id[2].strip()}"

    # 1) per-chain global alignment: crystal residue -> MD residue
    aligned, report, mdres = {}, [], {}
    for role, ch in roles.items():
        cres = [r for r in model[ch] if r.id[0] == " "]
        mres = md_residues(u, segids[ch])
        mdres[segids[ch]] = mres
        cs = "".join(seq1(r.get_resname()) for r in cres)
        ms = "".join(seq1(HIS.get(r.resname, r.resname)) for r in mres)
        pairs = align(cs, ms)
        mism = [(i, j) for i, j in pairs if i is not None and j is not None and cs[i] != ms[j]]
        crys_gap = [i for i, j in pairs if j is None]
        md_gap = [j for i, j in pairs if i is None]
        for i, j in pairs:
            if i is not None and j is not None:
                aligned[(ch, cres[i].id)] = (cres[i], mres[j])
        report.append(dict(pdb=pdb, role=role, chain=ch, segid=segids[ch], n_crystal=len(cres),
                           n_md=len(mres), n_aligned=len(pairs) - len(crys_gap) - len(md_gap),
                           mismatches=len(mism),
                           crystal_only=";".join(fmt(cres[i]) for i in crys_gap),
                           md_only=";".join(f"{mres[j].resname}{mres[j].resid}" for j in md_gap)))

    # 2) one rigid superposition of all aligned CA pairs, then per-residue coordinate verification
    cca, mca = [], []
    for c, m in aligned.values():
        if "CA" in c and "CA" in m.names:
            cca.append(c["CA"].coord); mca.append(X[m.index_of("CA")])
    R, t = kabsch(np.array(mca, float), np.array(cca, float))
    Xf = X @ R.T + t
    fit_rmsd = float(np.sqrt(((np.array(mca) @ R.T + t - np.array(cca)) ** 2).sum(1).mean()))

    def ca_dev(c, m):
        return float(np.linalg.norm(Xf[m.index_of("CA")] - c["CA"].coord))

    rows = []
    for kind, labs, slots in (("tcr", tlabs, tcr), ("mhc", mlabs, mhc)):
        for k, (lab, r) in enumerate(zip(labs, slots)):
            row = dict(pdb=pdb, kind=kind, slot=k, slot_label=lab)
            if r is None:
                row.update(status="absent_in_crystal")
                rows.append(row); continue
            ch = r.get_parent().id
            cheavy = [a for a in r if a.element != "H"]
            row.update(crystal_chain=ch, crystal_resnum=r.id[1], crystal_icode=r.id[2].strip(),
                       crystal_resname=r.get_resname(), n_heavy_crystal=len(cheavy))
            m = aligned.get((ch, r.id), (None, None))[1]
            note = ""
            if m is None or ca_dev(r, m) > COORD_TOL:
                cand = [x for x in mdres[segids[ch]] if "CA" in x.names
                        and seq1(HIS.get(x.resname, x.resname)) == seq1(r.get_resname())
                        and ca_dev(r, x) <= COORD_TOL]
                note = (f"alignment gave {m.resname}{m.resid}#{mdres[segids[ch]].index(m)}"
                        if m is not None else "alignment gap") + f"; coordinate match {len(cand)}"
                m = cand[0] if len(cand) == 1 else None
            if m is None:
                row.update(status="NOT_IN_MD", note=note)
                rows.append(row); continue
            mh = m.heavy
            d = np.linalg.norm(Xf[mh][None, :, :] - np.array([a.coord for a in cheavy])[:, None, :], axis=2)
            row.update(md_segid=m.segid, md_resid=m.resid, md_ordinal=mdres[m.segid].index(m),
                       md_resname=m.resname, n_heavy_md=len(mh),
                       ca_dev=round(ca_dev(r, m), 4), max_heavy_dev=round(float(d.min(1).max()), 4),
                       md_heavy_idx=",".join(map(str, mh)), note=note,
                       status="mapped" if seq1(r.get_resname()) == seq1(HIS.get(m.resname, m.resname))
                       else "RESNAME_MISMATCH")
            rows.append(row)
    return pd.DataFrame(rows), report, segids, roles, mdres, fit_rmsd


def groups_from_map_p84v2(df):
    """residue_map_<pdb>_p84v2.tsv rows -> (tcr_groups, mhc_groups) of MD heavy-atom indices."""
    def grp(kind, n):
        g = [None] * n
        for _, r in df[(df.kind == kind) & (df.status == "mapped")].iterrows():
            g[int(r.slot)] = np.array([int(x) for x in str(r.md_heavy_idx).split(",")])
        return g
    return grp("tcr", N_ROW), grp("mhc", N_COL)


def crystal_distance_matrix(pdb):
    """75 x 84 matrix of CRYSTAL-FRAME minimum heavy-atom distances (NaN where either residue is absent).

    Used for the distance-stratified rho of decision 4. Same residues, atom selection and function as the
    contact labels, so a cell at or below 4.5 A is exactly a crystal contact."""
    from contacts import PairMinDist
    from mdmap import crystal_model, chain_roles
    model = crystal_model(pdb)
    tcr, mhc = crystal_slots_p84v2(pdb, model, chain_roles(pdb))
    coords = []

    def groups(slots):
        g = []
        for r in slots:
            if r is None:
                g.append(None); continue
            idx = []
            for a in r:
                if a.element != "H":
                    idx.append(len(coords)); coords.append(a.coord.astype(float))
            g.append(np.array(idx))
        return g

    tg, mg = groups(tcr), groups(mhc)
    return PairMinDist(tg, mg)(np.array(coords))
