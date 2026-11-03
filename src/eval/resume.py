"""Generic resume-safe unit-of-work helper for the p82v1 rerun (and any future rerun).

Pattern: compute a unit's result, write it to a temp path, validate it, then atomically
rename (os.replace) to the unit's final cache path. A unit whose final cache path already
exists and passes validation is loaded and returned without recomputation. This gives
crash-safe, resumable execution without a separate job scheduler.
"""
import os
from pathlib import Path
import numpy as np


def _default_validate(d):
    for v in d.values():
        a = np.asarray(v)
        if a.dtype.kind in "fc" and not np.isfinite(a).all():
            return False
    return True


def resumable_unit(cache_dir, unit_id, compute_fn, validate_fn=None, verbose=True):
    """
    cache_dir : str or Path, created if missing.
    unit_id   : str, used as the cache filename stem (no extension, no path separators).
    compute_fn: () -> dict[str, array-like] to be cached.
    validate_fn: dict -> bool; defaults to "every float/complex array is all-finite".
    Returns (result_dict, was_resumed: bool).
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    final = cache_dir / f"{unit_id}.npz"
    vfn = validate_fn or _default_validate

    if final.exists():
        d = dict(np.load(final, allow_pickle=True))
        if vfn(d):
            if verbose: print(f"  [resume] {unit_id}: valid cache found, skipping recomputation")
            return d, True
        if verbose: print(f"  [resume] {unit_id}: cached file failed validation, recomputing")

    result = compute_fn()
    tmp = cache_dir / f"{unit_id}.tmp.npz"          # must end in .npz -- np.savez auto-appends
    np.savez(tmp, **result)                          # it otherwise, silently breaking the path
    d = dict(np.load(tmp, allow_pickle=True))
    if not vfn(d):
        raise RuntimeError(f"{unit_id}: freshly computed result failed validation, refusing to cache")
    os.replace(tmp, final)                            # atomic on POSIX filesystems
    if verbose: print(f"  [resume] {unit_id}: computed fresh, validated, cached at {final}")
    return d, False


def refuse_overwrite(path):
    """Call before writing any finalized, already-approved Stage A artifact."""
    if Path(path).exists():
        raise FileExistsError(
            f"refusing to overwrite existing finalized artifact {path} -- "
            f"bump the version tag (e.g. p82v2) for a corrected rerun instead")
