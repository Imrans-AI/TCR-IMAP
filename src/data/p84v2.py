"""p84v2 data layout: residues keyed by (IMGT number, insertion code), tokens placed by position.

New module (22 Oct 2026, results/review/CORRECTION_PLAN_P84V2.md). Nothing in src/models/grip.py or
any p82v1 caller is changed: `encode_seq` keeps its packed behaviour for existing callers.

Layout (from data/processed/pair_maps_p84v2.npz, built by scripts/53_build_pair_labels_p84v2.py):
  TCR axis  75 rows, one per (segment, IMGT number, insertion code) occupied in >= 1 complex, IMGT order
  MHC axis  84 columns, IMGT 50-90 and 1050-1090 plus 1061A and 1072A
Token rule (the alignment fix): the token at row r is the residue the label build assigns to row r
in that complex (`tcr_aa[c, r]`), 0 (padding) where the complex has no residue at that position.
"""
from pathlib import Path

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i + 1 for i, a in enumerate(AA)}                  # identical to src.models.grip
VOCAB = len(AA) + 1
SEG_ORDER = ["cdr1a", "cdr2a", "cdr3a", "cdr1b", "cdr2b", "cdr3b"]
GERMLINE_SEGS = {"cdr1a", "cdr2a", "cdr1b", "cdr2b"}
SLOT_START = {"cdr1": 27, "cdr2": 56, "cdr3": 105}
SLOT_LENS_P82V1 = {"cdr1a": 12, "cdr2a": 10, "cdr3a": 20, "cdr1b": 12, "cdr2b": 10, "cdr3b": 20}
POS82 = list(range(50, 91)) + list(range(1050, 1091))
ROOT = Path(__file__).resolve().parents[2]
NPZ = {"4p5A": ROOT / "data/processed/pair_maps_p84v2.npz", "6A": ROOT / "data/processed/pair_maps_p84v2_6A.npz"}


def imgt_rank(number, icode):
    """IMGT unique-numbering order. CDR3 insertions: 111 < 111A < 111B < ... < 112B < 112A < 112
    (111.x ascending, 112.x descending); elsewhere an insertion follows its base (1061 < 1061A)."""
    if not icode:
        return (int(number), 0, 0)
    k = ord(icode) - ord("A") + 1
    return (int(number), -1, -k) if int(number) == 112 else (int(number), 1, k)


def place_tokens(aa_rows):
    """THE token function of p84v2: one token per row, placed at that row (never packed).
    aa_rows: sequence of one-letter codes or '' (absent) per row -> int64 array, 0 = padding."""
    return np.array([AA_IDX.get(a, 0) if a else 0 for a in aa_rows], dtype=np.int64)


def place_tokens_batch(aa_matrix):
    return np.stack([place_tokens(r) for r in aa_matrix])


def load(cutoff="4p5A"):
    """Read the p84v2 label build. Returns dict with Y, S (float32), pdb, T (complex x 75), M
    (complex x 84), row/column manifests and the segment layout."""
    z = np.load(NPZ[cutoff], allow_pickle=True)
    T = place_tokens_batch(z["tcr_aa"]); M = place_tokens_batch(z["mhc_aa"])
    seg = [str(s) for s in z["row_segment"]]
    lens = {s: seg.count(s) for s in SEG_ORDER}
    assert seg == [s for s in SEG_ORDER for _ in range(lens[s])], "rows are not grouped by segment"
    return dict(Y=z["Y"].astype(np.float32), S=z["S"].astype(np.float32), pdb=[str(p) for p in z["pdb"]],
                T=T, M=M, seg_lens=lens, tcr_len=len(seg), mhc_len=M.shape[1],
                row_label=[str(x) for x in z["row_label"]], row_imgt_number=z["row_imgt_number"].astype(int),
                row_insertion_code=[str(x) for x in z["row_insertion_code"]],
                row_imgt_region=[str(x) for x in z["row_imgt_region"]],
                col_label=[str(x) for x in z["col_label"]])


_CACHE = {}


def _loaded(cutoff="4p5A"):
    if cutoff not in _CACHE:
        _CACHE[cutoff] = load(cutoff)
    return _CACHE[cutoff]


def T_for(pdbs, cutoff="4p5A", cdr3_last_imgt=None):
    """TCR token tensor (len(pdbs) x 75) for the given complexes, rows placed by position.
    cdr3_last_imgt: if given (AngleNet keeps CDR3 = IMGT 105-117), tokens of CDR3-segment rows with
    IMGT number above it are set to padding. Explicit, never applied by default."""
    d = _loaded(cutoff); pos = {p: i for i, p in enumerate(d["pdb"])}
    T = d["T"][[pos[str(p)] for p in pdbs]].copy()
    if cdr3_last_imgt is not None:
        drop = [r for r, (lab, num) in enumerate(zip(d["row_label"], d["row_imgt_number"]))
                if lab.startswith("cdr3") and num > cdr3_last_imgt]
        T[:, drop] = 0
    return T


def M_for(pdbs, cutoff="4p5A"):
    """MHC token tensor (len(pdbs) x 84), one token per column incl. 1061A and 1072A."""
    d = _loaded(cutoff); pos = {p: i for i, p in enumerate(d["pdb"])}
    return d["M"][[pos[str(p)] for p in pdbs]].copy()


def install_grip_shim(cutoff="4p5A"):
    """Process-local: make `src.models.grip` resolve to a copy whose layout constants are p84v2's
    (SEG_LENS 7/8/22/6/7/25, TCR_LEN 75), so model classes defined by exec'd p82v1 sources (19, 21,
    29b, 34b, 36b, 16b) are built for the p84v2 axes. Files are not changed; other processes and
    p82v1 callers are unaffected. encode_seq keeps its behaviour (peptides are still packed: they
    have no IMGT positions); TCR/MHC token construction is replaced in each runner and audited."""
    import sys
    import types
    import src.models.grip as real
    import src.models as pkg
    shim = types.ModuleType("src.models.grip")
    shim.__dict__.update({k: v for k, v in real.__dict__.items() if not k.startswith("__")})
    d = _loaded(cutoff)
    shim.SEG_LENS = dict(d["seg_lens"]); shim.TCR_LEN = d["tcr_len"]
    shim.SEG_ORDER = list(SEG_ORDER); shim.P84V2 = True
    shim.__file__ = real.__file__ + " [p84v2 shim]"
    sys.modules["src.models.grip"] = shim
    pkg.grip = shim
    return shim
