"""The one contact function of the MD occupancy analysis.

PairMinDist(tcr_groups, mhc_groups) holds, for each of the 84 TCR slots and 82 MHC positions, the
indices of that residue's heavy atoms in a coordinate array (None = slot not mapped). Calling it on
one frame's coordinates returns the 84 x 82 matrix of minimum heavy-atom distances (NaN where either
side is unmapped). Distances are minimum-image when a box is given (MDAnalysis distance_array with
the frame's box), so molecules split across the periodic boundary in the protein-only XTC need no
whole/nojump preprocessing; with box=None (crystal) they are plain Euclidean.

Contact = minimum heavy-atom distance <= CUTOFF (4.5 A), the p82v1 label definition (18b uses
`np.linalg.norm(a - b) <= CUTOFF` over all heavy-atom pairs, element != "H").
"""
import numpy as np
from MDAnalysis.lib.distances import distance_array

CUTOFF = 4.5


class PairMinDist:
    def __init__(self, tcr_groups, mhc_groups):
        self.nt, self.nm = len(tcr_groups), len(mhc_groups)
        self.t_slots = [k for k, g in enumerate(tcr_groups) if g is not None and len(g)]
        self.m_slots = [k for k, g in enumerate(mhc_groups) if g is not None and len(g)]
        self.a_idx = np.concatenate([np.asarray(tcr_groups[k]) for k in self.t_slots])
        self.b_idx = np.concatenate([np.asarray(mhc_groups[k]) for k in self.m_slots])
        # reduceat start offsets: atoms are concatenated slot by slot
        self.a_start = np.cumsum([0] + [len(tcr_groups[k]) for k in self.t_slots[:-1]])
        self.b_start = np.cumsum([0] + [len(mhc_groups[k]) for k in self.m_slots[:-1]])
        self._buf = np.empty((len(self.a_idx), len(self.b_idx)), dtype=np.float64)

    def __call__(self, X, box=None):
        """X: (n_atoms, 3) coordinates in A; box: MDAnalysis dimensions [lx, ly, lz, a, b, g] or None."""
        d = distance_array(X[self.a_idx], X[self.b_idx], box=box, result=self._buf)
        d = np.minimum.reduceat(np.minimum.reduceat(d, self.a_start, axis=0), self.b_start, axis=1)
        out = np.full((self.nt, self.nm), np.nan)
        out[np.ix_(self.t_slots, self.m_slots)] = d
        return out

    def contacts(self, X, box=None):
        """Boolean contact matrix (False where unmapped)."""
        d = self(X, box)
        return np.nan_to_num(d, nan=np.inf) <= CUTOFF


def groups_from_map(df):
    """residue_map_<pdb>.tsv rows -> (tcr_groups, mhc_groups) of MD heavy-atom indices."""
    def grp(kind, n):
        g = [None] * n
        for _, r in df[(df.kind == kind) & (df.status == "mapped")].iterrows():
            g[int(r.slot)] = np.array([int(x) for x in str(r.md_heavy_idx).split(",")])
        return g
    return grp("tcr", 84), grp("mhc", 82)
