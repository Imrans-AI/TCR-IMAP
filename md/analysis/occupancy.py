"""One pass over one trajectory: per-pair contact occupancy (all window frames, and docked window
frames) plus the per-frame TCR-pMHC contact count used by the docked-frame rule.

Definitions (md/MD_SYSTEM_SELECTION_CRITERIA.md, "Analysis definitions" and the 13 Oct secondary):
  * occupancy of a pair = fraction of window frames with minimum heavy-atom distance <= 4.5 A
    (contacts.PairMinDist, minimum image);
  * TCR-pMHC contact count = number of heavy-atom pairs within 4.5 A between TCR alpha + beta and
    peptide + MHC heavy chain (beta2m excluded), summed - the md/qc_replicate.sh definition;
  * docked frame = count >= DOCKED_FRAC x the count in the trajectory's first frame (t = 0).
Holds no scores and computes no correlation; scores enter only in stats.py at run time.
"""
import warnings

import numpy as np
from MDAnalysis.lib.distances import capped_distance

from contacts import PairMinDist, CUTOFF, groups_from_map
from mdmap import md_protein_universe

DOCKED_FRAC = 0.50


def is_docked(count, ref_count):
    """Docked frame: TCR-pMHC heavy-atom contact count >= 50 % of the t = 0 count (boundary included)."""
    return count >= DOCKED_FRAC * ref_count


def in_window(t_ps, t_start_ps, t_end_ps):
    """Analysis window, inclusive at both ends (20,000-200,000 ps -> 9,001 frames at 20 ps)."""
    return t_start_ps - 1e-3 <= t_ps <= t_end_ps + 1e-3


def heavy_indices(u, segids):
    ag = u.select_atoms("segid " + " ".join(segids))
    return ag.indices[~np.char.startswith(ag.names.astype(str), "H")]


def tcr_pmhc_pairs(X, box, tcr_idx, pmhc_idx):
    """Encoded (int64) heavy-atom contact pairs between TCR and peptide+MHC, i * n_pmhc + j."""
    p = capped_distance(X[tcr_idx], X[pmhc_idx], max_cutoff=CUTOFF, box=box, return_distances=False)
    return p[:, 0].astype(np.int64) * len(pmhc_idx) + p[:, 1].astype(np.int64)


def tcr_pmhc_count(X, box, tcr_idx, pmhc_idx):
    return len(tcr_pmhc_pairs(X, box, tcr_idx, pmhc_idx))


def is_docked_identity(persisting, n_ref_pairs):
    """Identity-based docked frame (criteria document, decision 3b, 24 Oct 2026): at least 50 % of the
    SPECIFIC heavy-atom contact pairs present at t = 0 are still present in the frame."""
    return persisting >= DOCKED_FRAC * n_ref_pairs


def run_trajectory(pdb, xtc, residue_map, chain_segments, t_start_ps, t_end_ps, docked=True,
                   groups_fn=groups_from_map, docked_identity=False):
    """Returns dict: occ (n_tcr x n_mhc, NaN unmapped), occ_docked (or None), n_window, n_docked,
    counts (per window frame), ref_count (t = 0), times (window), first_time.

    `groups_fn` maps a residue-map DataFrame to (tcr_groups, mhc_groups). It defaults to the p82v1
    84 x 82 reader, so existing callers are unchanged; the p84v2 analysis passes
    mdmap_p84v2.groups_from_map_p84v2 for the 75 x 84 axes. Nothing else in this function depends on
    the axis sizes.

    `docked_identity` additionally accumulates the identity-based docked-frame reading (decision 3b):
    the same per-frame contact-pair set is reused, so it costs nothing extra. It defaults to False so
    that existing callers behave exactly as before; the p84v2 production run passes True."""
    warnings.filterwarnings("ignore")
    u = md_protein_universe(pdb)
    u.load_new(str(xtc))
    f = PairMinDist(*groups_fn(residue_map))
    seg = chain_segments[chain_segments.pdb == pdb].set_index("role").segid
    tcr_idx = heavy_indices(u, [seg["TCRa"], seg["TCRb"]])
    pmhc_idx = heavy_indices(u, [seg["pep"], seg["MHC"]])

    acc = np.zeros((f.nt, f.nm)); acc_d = np.zeros((f.nt, f.nm)); acc_di = np.zeros((f.nt, f.nm))
    n = n_d = n_di = 0
    counts, persisting, times, ref, ref_pairs = [], [], [], None, None
    first_time = None
    for ts in u.trajectory:
        X = ts.positions.astype(np.float64)
        box = ts.dimensions
        if first_time is None:
            first_time = ts.time
            if docked:
                ref_pairs = tcr_pmhc_pairs(X, box, tcr_idx, pmhc_idx)
                ref = len(ref_pairs)
        if not in_window(ts.time, t_start_ps, t_end_ps):
            continue
        c = f.contacts(X, box)
        acc += c; n += 1
        times.append(ts.time)
        if docked:
            pk = tcr_pmhc_pairs(X, box, tcr_idx, pmhc_idx)
            k = len(pk)
            counts.append(k)
            if is_docked(k, ref):
                acc_d += c; n_d += 1
            if docked_identity:
                kept = int(np.isin(ref_pairs, pk).sum())
                persisting.append(kept)
                if is_docked_identity(kept, ref):
                    acc_di += c; n_di += 1
    unmapped = np.isnan(f(np.zeros((u.atoms.n_atoms, 3))))   # NaN pattern = unmapped slots
    occ = np.where(unmapped, np.nan, acc / max(n, 1))
    occ_d = np.where(unmapped, np.nan, acc_d / n_d) if (docked and n_d) else None
    occ_di = np.where(unmapped, np.nan, acc_di / n_di) if (docked_identity and n_di) else None
    return dict(occ=occ, occ_docked=occ_d, n_window=n, n_docked=n_d, counts=np.array(counts),
                ref_count=ref, times=np.array(times), first_time=first_time,
                occ_docked_identity=occ_di, n_docked_identity=n_di,
                persisting=np.array(persisting), n_ref_pairs=ref)
