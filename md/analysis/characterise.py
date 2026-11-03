"""Simulation-quality and interface-characterisation measures, fixed in
md/MD_SYSTEM_SELECTION_CRITERIA.md ("Analysis decisions settled, and the characterisation measures added -
24 October 2026").

**These are score-independent.** Nothing in this module reads, imports or accepts a model prediction. It is
called from phase A of the pipeline, before any score file is opened, which is what keeps occupancy and the
characterisation measures fixed before any score exists.

Measures, all over the 20-200 ns window unless stated:
  * backbone_rmsd       - protein Ca RMSD against the energy-minimised starting structure (md/runs/<pdb>/em/em.gro)
  * interface_rmsd      - Ca RMSD of the TCR CDR-loop residues after a least-squares fit on the MHC alpha1/alpha2
                          Ca atoms, relative to the crystal frame (step3_input, which carries the crystal pose)
  * per-loop RMSF       - per-residue Ca RMSF after the same fit, averaged within each of the six CDR loops
  * buried_sasa         - SASA(TCR) + SASA(pMHC) - SASA(complex), Shrake-Rupley, probe 1.4 A
  * hydrogen bonds      - donor-acceptor heavy-atom distance <= 3.5 A and D-H...A angle >= 150 deg, counted
                          TCR<->MHC alpha1/alpha2 and TCR<->peptide separately

CDR residues for the loop measures are the canonical loops - IMGT 27-38 (CDR1), 56-65 (CDR2) and **105-117**
(CDR3) - so the seven J-region framework rows per chain (IMGT 118-124) that the pair universe carries are
excluded from RMSF and interface RMSD. They are framework, not CDR. The pair universe itself is unchanged.

SASA implementation: Shrake-Rupley with a golden-spiral sphere of SASA_N_POINTS points and probe radius
SASA_PROBE, over heavy atoms, with element radii VDW_RADII. Buried SASA is computed exactly but only over the
atoms that can contribute: an atom's SASA changes between "its own group alone" and "the whole complex" only if
some atom of the other group lies within r_i + r_j + 2 * probe of it, so the sum is restricted to those atoms.
`buried_sasa_naive` computes the same quantity by the definition and is used in the tests to check this.

Periodic boundaries (24 October 2026). Every coordinate-based measure here needs a trajectory made whole and
jump-free (md/analysis/pbc.py, the procedure of md/qc_replicate.sh). characterise_trajectory checks every window frame
with `assembly_problem` and raises on a broken or unassembled complex rather than return a wrong number.
"""
import numpy as np
from scipy.spatial import cKDTree

SASA_PROBE = 1.4          # A, standard water probe
SASA_N_POINTS = 256       # golden-spiral points per atom
VDW_RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "H": 1.20}
VDW_DEFAULT = 1.70
HB_D_A_CUTOFF = 3.5       # A, donor-acceptor heavy-atom distance
HB_ANGLE_CUTOFF = 150.0   # degrees, D-H...A
CDR_RANGES = {"cdr1": (27, 38), "cdr2": (56, 65), "cdr3": (105, 117)}   # canonical loops; FR4 118-124 excluded


# ----------------------------------------------------------------- geometry helpers
def sphere_points(n=SASA_N_POINTS):
    """n roughly equidistant unit vectors (golden spiral)."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])


def radii_for(names, elements=None):
    """Per-atom van der Waals radius from the element (first letter of the atom name when no element)."""
    out = []
    for k, nm in enumerate(names):
        e = (elements[k] if elements is not None and str(elements[k]).strip() else str(nm).strip()[:1]).upper()
        out.append(VDW_RADII.get(e, VDW_DEFAULT))
    return np.asarray(out, float)


def _sasa_of(X, R, which, tree, sphere):
    """Per-atom SASA of the atoms `which`, in the environment described by (X, R, tree)."""
    out = np.zeros(len(which))
    rmax = R.max()
    for n, i in enumerate(which):
        Ri = R[i] + SASA_PROBE
        nb = np.array(tree.query_ball_point(X[i], Ri + rmax + SASA_PROBE), dtype=int)
        nb = nb[nb != i]
        pts = X[i] + Ri * sphere
        if len(nb):
            d2 = ((pts[:, None, :] - X[nb][None, :, :]) ** 2).sum(-1)
            acc = (d2 >= ((R[nb] + SASA_PROBE) ** 2)[None, :]).all(axis=1)
        else:
            acc = np.ones(len(sphere), bool)
        out[n] = 4 * np.pi * Ri ** 2 * acc.mean()
    return out


def buried_sasa(X, R, idx_a, idx_b, sphere=None):
    """SASA(a) + SASA(b) - SASA(a+b), computed exactly.

    Only atoms with an atom of the other group within r_i + r_j + 2*probe can differ between the isolated
    group and the complex, so the difference is summed over those atoms alone."""
    sphere = sphere_points() if sphere is None else sphere
    Xa, Xb = X[idx_a], X[idx_b]
    ta, tb = cKDTree(Xa), cKDTree(Xb)
    reach = R.max() + R.max() + 2 * SASA_PROBE
    pairs = ta.query_ball_tree(tb, reach)
    near_a = np.array([i for i, pl in enumerate(pairs) if pl], dtype=int)
    near_b = np.array(sorted({j for pl in pairs for j in pl}), dtype=int)
    if not len(near_a):
        return 0.0
    Xab = np.concatenate([Xa, Xb]); Rab = np.concatenate([R[idx_a], R[idx_b]])
    tab = cKDTree(Xab)
    sa_alone = _sasa_of(Xa, R[idx_a], near_a, ta, sphere)
    sb_alone = _sasa_of(Xb, R[idx_b], near_b, tb, sphere)
    sa_cplx = _sasa_of(Xab, Rab, near_a, tab, sphere)
    sb_cplx = _sasa_of(Xab, Rab, near_b + len(Xa), tab, sphere)
    return float((sa_alone - sa_cplx).sum() + (sb_alone - sb_cplx).sum())


def buried_sasa_naive(X, R, idx_a, idx_b, sphere=None):
    """The same quantity straight from the definition (every atom, three full SASA computations). Tests only."""
    sphere = sphere_points() if sphere is None else sphere
    def total(idx):
        Xi, Ri = X[idx], R[idx]
        return _sasa_of(Xi, Ri, np.arange(len(idx)), cKDTree(Xi), sphere).sum()
    both = np.concatenate([idx_a, idx_b])
    return float(total(idx_a) + total(idx_b) - total(both))


def kabsch_fit(P, Q):
    """R, t minimising |P @ R.T + t - Q| (rows = points)."""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    Rm = Vt.T @ D @ U.T
    return Rm, qc - pc @ Rm.T


def rmsd(A, B):
    return float(np.sqrt(((A - B) ** 2).sum(1).mean()))


def fitted_rmsd(X, ref, fit_idx, measure_idx):
    """Superpose X on ref using fit_idx, then RMSD over measure_idx."""
    Rm, t = kabsch_fit(X[fit_idx], ref[fit_idx])
    return rmsd(X[measure_idx] @ Rm.T + t, ref[measure_idx])


class RMSFAccumulator:
    """Mean and mean-square of fitted Ca positions, in one pass (RMSF = sqrt(<|x|^2> - |<x>|^2))."""
    def __init__(self, n):
        self.n = 0
        self.s = np.zeros((n, 3))
        self.s2 = np.zeros(n)

    def add(self, Xfit):
        self.s += Xfit
        self.s2 += (Xfit ** 2).sum(1)
        self.n += 1

    def rmsf(self):
        mean = self.s / self.n
        return np.sqrt(np.maximum(self.s2 / self.n - (mean ** 2).sum(1), 0.0))


# ----------------------------------------------------------------- periodic-boundary guard
MAX_CA_STEP = 4.5         # A: consecutive Ca atoms within one chain (3.8 trans, 2.9 cis); larger = chain broken
MAX_CHAIN_GAP = 15.0      # A: each chain's nearest Ca to the rest of the complex; the five chains must be connected


def assembly_problem(X, chain_ca):
    """None if the frame's coordinates hold every chain whole and the complex assembled, without periodic wrapping;
    otherwise a description. These measures use raw coordinates, so a chain broken across the box or sitting in
    another periodic image corrupts them (criteria document, 24 October 2026); this makes that fail loudly."""
    for seg, idx in chain_ca.items():
        steps = np.linalg.norm(np.diff(X[idx], axis=0), axis=1)
        if len(steps) and steps.max() > MAX_CA_STEP:
            return f"chain {seg} broken: consecutive Ca {steps.max():.1f} A apart"
    segs = list(chain_ca)
    trees = {s: cKDTree(X[chain_ca[s]]) for s in segs}
    linked = {segs[0]}
    grew = True
    while grew:
        grew = False
        for s in segs:
            if s in linked:
                continue
            if any(trees[s].sparse_distance_matrix(trees[t], MAX_CHAIN_GAP).nnz for t in linked):
                linked.add(s); grew = True
    if len(linked) < len(segs):
        return f"complex not assembled: {sorted(set(segs) - linked)} more than {MAX_CHAIN_GAP} A from the rest"
    return None


# ----------------------------------------------------------------- MD selections and the trajectory pass
def selections(pdb, residue_map, chain_segments, u):
    """Atom-index selections for one system, from the p84v2 residue map and the chain segments.

    Returns dict with: all_ca (protein), mhc_ca (the mapped alpha1/alpha2 columns), loop_ca {loop: idx},
    loop_resid {loop: labels}, tcr_heavy, pmhc_heavy, pep_heavy, mhc12_heavy."""
    import re
    names = u.atoms.names.astype(str)
    seg = chain_segments[chain_segments.pdb == pdb].set_index("role").segid

    def heavy(segids):
        ag = u.select_atoms("segid " + " ".join(segids))
        return ag.indices[~np.char.startswith(ag.names.astype(str), "H")]

    def ca_of(row):
        idx = np.array([int(x) for x in str(row.md_heavy_idx).split(",")])
        ca = idx[names[idx] == "CA"]
        return int(ca[0]) if len(ca) == 1 else None

    mapped = residue_map[residue_map.status == "mapped"]
    mhc_rows = mapped[mapped.kind == "mhc"]
    mhc_ca = np.array([c for c in (ca_of(r) for r in mhc_rows.itertuples()) if c is not None])

    loop_ca, loop_lab = {}, {}
    for r in mapped[mapped.kind == "tcr"].itertuples():
        m = re.match(r"(cdr[123])([ab]):(\d+)([A-Z]?)$", str(r.slot_label))
        assert m, r.slot_label
        loop, chain, num = m.group(1), m.group(2), int(m.group(3))
        lo, hi = CDR_RANGES[loop]
        if not (lo <= num <= hi):          # J-region framework rows (IMGT 118-124) are not CDR
            continue
        c = ca_of(r)
        if c is None:
            continue
        key = f"{loop}{chain}"
        loop_ca.setdefault(key, []).append(c)
        loop_lab.setdefault(key, []).append(str(r.slot_label))
    loop_ca = {k: np.array(v) for k, v in sorted(loop_ca.items())}

    return dict(all_ca=u.select_atoms("name CA").indices,
                mhc_ca=mhc_ca, loop_ca=loop_ca, loop_resid=loop_lab,
                tcr_heavy=heavy([seg["TCRa"], seg["TCRb"]]),
                pmhc_heavy=heavy([seg["pep"], seg["MHC"]]),
                pep_heavy=heavy([seg["pep"]]),
                mhc12_heavy=np.concatenate([np.array([int(x) for x in str(r.md_heavy_idx).split(",")])
                                            for r in mhc_rows.itertuples()]))


def em_reference_ca(pdb, u, root):
    """Protein Ca coordinates of the energy-minimised starting structure (md/runs/<pdb>/em/em.gro),
    in the universe's atom order."""
    import MDAnalysis as mda
    g = mda.Universe(str(root / f"md/runs/{pdb}/em/em.gro"))
    prot = g.select_atoms("protein")
    assert prot.n_atoms >= u.atoms.n_atoms, (prot.n_atoms, u.atoms.n_atoms)
    sub = g.atoms[:u.atoms.n_atoms]
    assert list(sub.names.astype(str)) == list(u.atoms.names.astype(str)), \
        "em.gro atom order differs from the protein-only universe"
    return sub.positions.copy()


def characterise_trajectory(pdb, xtc, residue_map, chain_segments, t_start_ps, t_end_ps, root,
                            sasa_stride=50, log=print):
    """One pass for RMSD / interface RMSD / RMSF, a strided pass for buried SASA. No scores.

    Returns dict with per-frame series and summaries."""
    import warnings
    warnings.filterwarnings("ignore")
    from mdmap import md_protein_universe
    from occupancy import in_window
    u = md_protein_universe(pdb)
    sel = selections(pdb, residue_map, chain_segments, u)
    crystal = u.atoms.positions.copy()                 # step3_input carries the crystal pose
    em_ca = em_reference_ca(pdb, u, root)
    u.load_new(str(xtc))
    chain_ca = {s: u.select_atoms(f"segid {s} and name CA").indices for s in sorted(set(u.atoms.segids))}

    loops = sorted(sel["loop_ca"])
    loop_idx = np.concatenate([sel["loop_ca"][k] for k in loops])
    loop_of = np.concatenate([[k] * len(sel["loop_ca"][k]) for k in loops])
    acc = RMSFAccumulator(len(loop_idx))
    sphere = sphere_points()
    names = u.atoms.names.astype(str)
    elements = getattr(u.atoms, "elements", None)
    R_all = radii_for(names, np.asarray(elements).astype(str) if elements is not None else None)

    times, bb, irmsd, sasa_t, sasa_v = [], [], [], [], []
    for n, ts in enumerate(u.trajectory):
        if not in_window(ts.time, t_start_ps, t_end_ps):
            continue
        X = ts.positions.astype(np.float64)
        bad = assembly_problem(X, chain_ca)
        if bad:
            raise ValueError(f"{pdb} {xtc}: t = {ts.time:.0f} ps: {bad}. Coordinate measures need a whole, "
                             f"jump-free trajectory (md/analysis/pbc.py)")
        bb.append(fitted_rmsd(X, em_ca, sel["all_ca"], sel["all_ca"]) if em_ca.shape == X.shape
                  else fitted_rmsd(X[sel["all_ca"]], em_ca[sel["all_ca"]],
                                   np.arange(len(sel["all_ca"])), np.arange(len(sel["all_ca"]))))
        irmsd.append(fitted_rmsd(X, crystal, sel["mhc_ca"], loop_idx))
        Rm, t = kabsch_fit(X[sel["mhc_ca"]], crystal[sel["mhc_ca"]])
        acc.add(X[loop_idx] @ Rm.T + t)
        times.append(ts.time)
        if len(times) % sasa_stride == 1 or sasa_stride == 1:
            sasa_t.append(ts.time)
            sasa_v.append(buried_sasa(X, R_all, sel["tcr_heavy"], sel["pmhc_heavy"], sphere))
    rmsf = acc.rmsf()
    per_loop = {k: float(rmsf[loop_of == k].mean()) for k in loops}
    log(f"  {pdb} {xtc.name}: {len(times)} window frames, SASA on {len(sasa_v)} (stride {sasa_stride})")
    return dict(times=np.array(times), backbone_rmsd=np.array(bb), interface_rmsd=np.array(irmsd),
                rmsf_per_residue=rmsf, rmsf_labels=np.concatenate([sel["loop_resid"][k] for k in loops]),
                rmsf_loop=loop_of, per_loop_rmsf=per_loop,
                sasa_times=np.array(sasa_t), buried_sasa=np.array(sasa_v), sasa_stride=sasa_stride,
                n_window=len(times))


def hbond_counts(pdb, xtc, residue_map, chain_segments, t_start_ps, t_end_ps, stride=1):
    """Interface hydrogen bonds: TCR<->MHC alpha1/alpha2 and TCR<->peptide, counted separately.
    Criteria: donor-acceptor heavy-atom distance <= 3.5 A, D-H...A angle >= 150 deg. No scores."""
    import warnings
    warnings.filterwarnings("ignore")
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBA
    from mdmap import md_protein_universe
    from occupancy import in_window
    u = md_protein_universe(pdb)
    sel = selections(pdb, residue_map, chain_segments, u)
    u.load_new(str(xtc))
    frames = [i for i, ts in enumerate(u.trajectory) if in_window(ts.time, t_start_ps, t_end_ps)]
    frames = frames[::stride]

    def idx_sel(idx):
        return "index " + " ".join(map(str, idx))

    out = {}
    for tag, other in (("TCR-MHC", sel["mhc12_heavy"]), ("TCR-peptide", sel["pep_heavy"])):
        h = HBA(universe=u, between=[idx_sel(sel["tcr_heavy"]), idx_sel(other)],
                d_a_cutoff=HB_D_A_CUTOFF, d_h_a_angle_cutoff=HB_ANGLE_CUTOFF, update_selections=False)
        h.run(frames=frames, verbose=False)
        res = h.results.hbonds
        counts = np.array([int((res[:, 0] == f).sum()) for f in np.unique(res[:, 0])]) if len(res) else np.zeros(0)
        per_frame = np.zeros(len(frames))
        if len(res):
            uf, cf = np.unique(res[:, 0], return_counts=True)
            pos = {int(f): k for k, f in enumerate(frames)}
            for f, c in zip(uf, cf):
                per_frame[pos[int(f)]] = c
        pair_key = {}
        if len(res):
            dres = u.atoms[res[:, 1].astype(int)].resids
            ares = u.atoms[res[:, 3].astype(int)].resids
            dseg = u.atoms[res[:, 1].astype(int)].segids
            aseg = u.atoms[res[:, 3].astype(int)].segids
            for ds, dr, as_, ar in zip(dseg, dres, aseg, ares):
                pair_key[f"{ds}{dr}-{as_}{ar}"] = pair_key.get(f"{ds}{dr}-{as_}{ar}", 0) + 1
        out[tag] = dict(n_frames=len(frames), mean_per_frame=float(per_frame.mean()),
                        sd_per_frame=float(per_frame.std()), per_frame=per_frame,
                        pair_occupancy={k: v / len(frames) for k, v in pair_key.items()})
    return out
