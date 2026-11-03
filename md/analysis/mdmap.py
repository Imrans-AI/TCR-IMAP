"""Residue mapping for the MD occupancy analysis: each of the 84 TCR slots and 82 MHC positions of the
p82v1 labels -> one residue of the CHARMM-GUI topology, by global sequence alignment, never by number.

Crystal side reproduces scripts/18b_build_pair_labels_p82v1.py exactly:
  * residues are keyed by IMGT number alone, standard residues only (`{x.id[1]: x for x in chain if
    x.id[0] == " "}`), so where an insertion code repeats a number the residue LAST in file order
    wins (e.g. MHC slot 1061 -> residue 1061A, slot 1072 -> 1072A in all four systems);
  * TCR slot k of segment `seg`+`tag` is IMGT position lo + k for k < SEG_LENS[key], with
    SEG_RANGE = cdr1 27, cdr2 56, cdr3 105 (so cdr3 slots run 105-124);
  * MHC columns are IMGT 50-90 and 1050-1090 on `mhc_chain1`.
MD side: the CHARMM-GUI segment built from each crystal chain (identified by the chain's own
`<pdb>_clean_pro<x>.pdb`, which carries the segid), read from gromacs/step3_input.psf/.pdb.
MD residues are NOT taken from MDAnalysis' residue grouping: CHARMM-GUI drops insertion codes, so
crystal 1061 and 1061A both become resid 1061 ALA and MDAnalysis merges them into one 20-atom residue.
`md_residues()` splits a segment into residues wherever resid or resname changes or an atom name
repeats. Every alignment-mapped residue is then verified by coordinates: step3_input keeps the
crystal coordinates (rigid offset only), so after one superposition each mapped MD residue's heavy
atoms must coincide with the crystal residue's (CA within COORD_TOL), which also resolves any
alignment ambiguity between identical neighbours (e.g. Ala1061/Ala1061A).

Reads no trajectory, no score and no label value other than S/Y for cross-checks elsewhere.
"""
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio import Align

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.models.grip import SEG_ORDER, SEG_LENS, TCR_LEN  # noqa: E402

SYSTEMS = ["5yxn", "7n1f", "3vxs", "8eo8"]          # queue order, H L H L
SEG_RANGE = {"cdr1": (27, 38), "cdr2": (56, 65), "cdr3": (105, 117)}   # as 18b
POS82 = list(range(50, 91)) + list(range(1050, 1091))                  # as 18b
CUTOFF = 4.5                                                           # as 18b / label definition
HIS = {"HSD": "HIS", "HSE": "HIS", "HSP": "HIS"}
PROTEIN_SEGS = ["PROA", "PROB", "PROC", "PROD", "PROE"]
COORD_TOL = 0.05   # A, CA position after superposition (input coordinates are crystal coordinates)


def gromacs_dir(pdb):
    return Path(glob.glob(str(ROOT / f"md/systems/{pdb}/charmm-gui-*/gromacs"))[0])


def chain_roles(pdb):
    """role -> crystal chain id, from the STCRDab index used by 18b."""
    idx = pd.read_csv(ROOT / "data/interim/stcrdab_index.tsv", sep="\t").drop_duplicates("pdb").set_index("pdb")
    r = idx.loc[pdb]
    return {"TCRa": str(r.Achain), "TCRb": str(r.Bchain), "MHC": str(r.mhc_chain1),
            "b2m": str(r.mhc_chain2), "pep": str(r.antigen_chain)}


def crystal_model(pdb):
    return PDBParser(QUIET=True).get_structure(pdb, str(ROOT / f"md/systems/{pdb}/prep/{pdb}_clean.pdb"))[0]


def slot_labels():
    """The 84 TCR slot labels in 18b order ('cdr1a:27' ...) and the 82 MHC positions."""
    labs = []
    for tag in ("a", "b"):
        for seg, (lo, _) in SEG_RANGE.items():
            key = f"{seg}{tag}"
            labs += [f"{key}:{p}" for p in range(lo, lo + SEG_LENS[key])]
    assert len(labs) == TCR_LEN
    return labs, list(POS82)


def crystal_slots(model, roles):
    """18b's slot -> residue logic, verbatim. Returns (tcr, mhc): lists of Bio residues or None."""
    tcr = [None] * TCR_LEN
    for tag, role in (("a", "TCRa"), ("b", "TCRb")):
        tres = {x.id[1]: x for x in model[roles[role]] if x.id[0] == " "}     # last-wins, as 18b
        for seg, (lo, _) in SEG_RANGE.items():
            key = f"{seg}{tag}"
            off = sum(SEG_LENS[k] for k in SEG_ORDER[:SEG_ORDER.index(key)])
            for p in range(lo, lo + SEG_LENS[key]):
                if p in tres:
                    tcr[off + (p - lo)] = tres[p]
    mres = {x.id[1]: x for x in model[roles["MHC"]] if x.id[0] == " "}       # last-wins, as 18b
    mhc = [mres.get(p) for p in POS82]
    return tcr, mhc


def md_protein_universe(pdb):
    """Protein-only MDAnalysis Universe in SOLU atom order (= the XTC atom order)."""
    import MDAnalysis as mda
    G = gromacs_dir(pdb)
    u = mda.Universe(str(G / "step3_input.psf"), str(G / "step3_input.pdb"))
    prot = u.select_atoms("segid " + " ".join(PROTEIN_SEGS))
    solu = [int(x) for x in open(G / "index.ndx").read().split("[")[1].split("]")[1].split()]
    assert list(prot.indices + 1) == solu, "protein atoms are not the SOLU group in order"
    return mda.Merge(prot)


def chain_segids(pdb, model):
    """crystal chain id -> CHARMM-GUI segid, from CHARMM-GUI's per-chain input files (exact
    residue-name sequence match), not from file order."""
    out = {}
    cg = gromacs_dir(pdb).parent
    for f in sorted(cg.glob(f"{pdb}_clean_pro?.pdb")):
        seq, segid = [], None
        last = None
        for line in open(f):
            if line.startswith("ATOM"):
                segid = line[72:76].strip()
                key = line[22:27]
                if key != last:
                    seq.append(HIS.get(line[17:21].strip(), line[17:21].strip()))
                    last = key
        hits = [c.id for c in model if [r.get_resname() for r in c if r.id[0] == " "] == seq]
        assert len(hits) == 1, f"{f.name}: matches crystal chains {hits}"
        out[hits[0]] = segid
    assert len(set(out.values())) == len(out)
    return out


class MDRes:
    """One MD residue: a contiguous block of atoms in SOLU order."""
    def __init__(self, segid, resid, resname, indices, names):
        self.segid, self.resid, self.resname = segid, int(resid), resname
        self.indices, self.names = np.asarray(indices), list(names)

    @property
    def heavy(self):
        return np.array([i for i, n in zip(self.indices, self.names) if not n.startswith("H")])

    def index_of(self, name):
        return int(self.indices[self.names.index(name)])


def md_residues(u, segid):
    """Split one segment into residues: new residue when resid or resname changes or an atom name
    repeats (two consecutive same-name residues sharing a resid, e.g. Ala1061/Ala1061A)."""
    ag = u.select_atoms(f"segid {segid}")
    out, cur, seen = [], None, set()
    for i, n, rid, rn in zip(ag.indices, ag.names, ag.resids, ag.resnames):
        if cur is None or rid != cur[1] or rn != cur[2] or n in seen:
            if cur is not None:
                out.append(MDRes(segid, cur[1], cur[2], cur[3], cur[4]))
            cur, seen = [segid, rid, rn, [], []], set()
        cur[3].append(i); cur[4].append(n); seen.add(n)
    out.append(MDRes(segid, cur[1], cur[2], cur[3], cur[4]))
    return out


def align(crys_seq, md_seq):
    """Global alignment; end gaps free, internal gaps penalised. Returns list of (i_crys, j_md)
    aligned index pairs (None on a gap side)."""
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.match_score, al.mismatch_score = 2, -3
    al.open_gap_score, al.extend_gap_score = -5, -1
    al.end_gap_score = 0
    a = al.align(crys_seq, md_seq)[0]
    pairs, i, j = [], 0, 0
    for (cs, ce), (ms, me) in zip(*a.aligned):
        while i < cs: pairs.append((i, None)); i += 1
        while j < ms: pairs.append((None, j)); j += 1
        for k in range(ce - cs):
            pairs.append((cs + k, ms + k))
        i, j = ce, me
    while i < len(crys_seq): pairs.append((i, None)); i += 1
    while j < len(md_seq): pairs.append((None, j)); j += 1
    return pairs


def heavy(names):
    return [n for n in names if not n.startswith("H")]


def kabsch(P, Q):
    """R, t minimising |P @ R.T + t - Q| (rows = points)."""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    return R, qc - pc @ R.T


def build_map(pdb):
    """Residue map + gate G1 facts for one system.
    Returns (map DataFrame, per-chain report list, chain->segid, role->chain, MD residues by segid)."""
    model = crystal_model(pdb)
    roles = chain_roles(pdb)
    segids = chain_segids(pdb, model)
    u = md_protein_universe(pdb)
    X = u.atoms.positions
    tcr, mhc = crystal_slots(model, roles)
    tlabs, mpos = slot_labels()
    fmt = lambda r: f"{r.get_resname()}{r.id[1]}{r.id[2].strip()}"

    # 1) per-chain global alignment: crystal residue (Bio) -> MD residue (MDRes)
    aligned, report, mdres = {}, [], {}
    for role, ch in roles.items():
        cres = [r for r in model[ch] if r.id[0] == " "]
        mres = md_residues(u, segids[ch])
        mdres[segids[ch]] = mres
        cs = "".join(seq1(r.get_resname()) for r in cres)
        ms = "".join(seq1(HIS.get(r.resname, r.resname)) for r in mres)
        pairs = align(cs, ms)
        mism = [(i, j) for i, j in pairs if i is not None and j is not None and cs[i] != ms[j]]
        crys_gap = [i for i, j in pairs if j is None]       # crystal residue absent from MD
        md_gap = [j for i, j in pairs if i is None]          # MD residue absent from crystal (modelled)
        for i, j in pairs:
            if i is not None and j is not None:
                aligned[(ch, cres[i].id)] = (cres[i], mres[j])
        report.append(dict(pdb=pdb, role=role, chain=ch, segid=segids[ch], n_crystal=len(cres),
                           n_md=len(mres), n_aligned=len(pairs) - len(crys_gap) - len(md_gap),
                           mismatches=len(mism),
                           crystal_only=";".join(fmt(cres[i]) for i in crys_gap),
                           md_only=";".join(f"{mres[j].resname}{mres[j].resid}" for j in md_gap)))

    # 2) coordinate verification: one superposition of all aligned CA pairs, then per residue
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
    for kind, labs, slots in (("tcr", tlabs, tcr), ("mhc", [str(p) for p in mpos], mhc)):
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
                # resolve by coordinates within the same segment (same residue type only)
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
            # heavy-atom coordinate check: every crystal heavy atom has an MD heavy atom within tol
            d = np.linalg.norm(Xf[mh][None, :, :] - np.array([a.coord for a in cheavy])[:, None, :], axis=2)
            row.update(md_segid=m.segid, md_resid=m.resid, md_ordinal=mdres[m.segid].index(m),
                       md_resname=m.resname, n_heavy_md=len(mh),
                       ca_dev=round(ca_dev(r, m), 4), max_heavy_dev=round(float(d.min(1).max()), 4),
                       md_heavy_idx=",".join(map(str, mh)), note=note,
                       status="mapped" if seq1(r.get_resname()) == seq1(HIS.get(m.resname, m.resname))
                       else "RESNAME_MISMATCH")
            rows.append(row)
    return pd.DataFrame(rows), report, segids, roles, mdres, fit_rmsd


if __name__ == "__main__":
    print(slot_labels()[0][:3])
