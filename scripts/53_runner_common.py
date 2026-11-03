"""Shared machinery for the p84v2 consumer runners (scripts/53_*_p84v2.py).

A runner executes an existing p82v1 consumer script's OWN source, unchanged in its models,
hyperparameters, folds, metrics and statistics, with three kinds of explicit, counted edits:
  1. TCR and MHC token tensors: the packed `encode_seq` constructions are replaced by
     src.data.p84v2.T_for / M_for (tokens placed at the row the manifest assigns each residue);
  2. data paths: pair_maps/pair_meta p82v1 -> p84v2, tcr_graphs.npz -> tcr_graphs_p84v2.npz;
  3. output paths: *_p82v1.* and _partial/p82v1/ -> p84v2 (nothing p82v1 can be written).
Layout constants reach every model class (including classes defined by nested exec of 19/21)
through src.data.p84v2.install_grip_shim(), process-local. After patching, the source is scanned and
the runner refuses to continue if any p82v1 data/output path or packed TCR/MHC token call remains.
Modes: --check (default) executes the patched definitions, calls the consumer's own data path up to
the token tensors (early-return hook) and instantiates its models - no training, no inference.
--run executes main(); it is NOT used by the rebuild task (CORRECTION_PLAN_P84V2.md: stop before training).
"""
import argparse
import hashlib
import inspect
import re
import sys
from pathlib import Path

sys.path.insert(0, ".")
from src.data import p84v2 as P  # noqa: E402

PATCHED_DIR = Path("results/review/rebuild_p84v2/patched_sources")
FORBIDDEN = [r"encode_seq\(r\[k\]", r"\.mhc_contact\]", r"pair_maps_p82v1", r"pair_meta_p82v1", r"tcr_graphs\.npz",
             r"_partial/p82v1", r"_p82v1\.(csv|npz|jsonl)", r"_p82v1_(4p5A|6A)"]
T_RE = re.compile(r"(?m)^(?P<ind>[ \t]*)T\s*=\s*np\.stack\(\[np\.concatenate\(\[encode_seq\(r\[k\],\s*SEG_LENS\[k\]\) "
                  r"for k in SEG_ORDER\]\) for _,r in (?P<frame>meta|df)\.iterrows\(\)\]\)")
M_RE = re.compile(r"(?m)^(?P<ind>[ \t]*)M\s*=\s*np\.stack\(\[encode_seq\(s,\s*[^\n]*?\) for s in (?P<frame>meta|df)\.mhc_contact\]\)")
PATH_REPS = [("pair_maps_p82v1", "pair_maps_p84v2"), ("pair_meta_p82v1", "pair_meta_p84v2"),
             ("data/processed/tcr_graphs.npz", "data/processed/tcr_graphs_p84v2.npz"),
             ("_partial/p82v1/", "_partial/p84v2/"), ("_p82v1.csv", "_p84v2.csv"), ("_p82v1.npz", "_p84v2.npz"),
             ("_p82v1.jsonl", "_p84v2.jsonl"), ("_p82v1_4p5A", "_p84v2_4p5A"), ("_p82v1_6A", "_p84v2_6A")]


def patch(src, extra=(), cdr3_last_imgt=None, hook_after=None):
    """Return (patched source, edit log). Every edit is counted; nothing is silent."""
    log = []
    arg = "" if cdr3_last_imgt is None else f", cdr3_last_imgt={cdr3_last_imgt}"
    src, n = T_RE.subn(lambda m: f"{m['ind']}T = _P84.T_for({m['frame']}.pdb{arg})", src); log.append(("TCR tokens -> _P84.T_for", n))
    src, n = M_RE.subn(lambda m: f"{m['ind']}M = _P84.M_for({m['frame']}.pdb)", src); log.append(("MHC tokens -> _P84.M_for", n))
    for old, new, want in extra:
        k = src.count(old)
        assert k == want, f"expected {want} occurrence(s) of {old!r}, found {k}"
        src = src.replace(old, new); log.append((f"{old!r} -> {new!r}", k))
    for old, new in PATH_REPS:
        k = src.count(old); src = src.replace(old, new)
        if k: log.append((f"{old} -> {new}", k))
    if hook_after is not None:
        rx = re.compile(hook_after, re.M); m = list(rx.finditer(src))
        assert len(m) == 1, f"hook anchor {hook_after!r} matched {len(m)} times"
        ind = re.match(r"[ \t]*", src[m[0].start():]).group(0)
        src = src[:m[0].end()] + f"\n{ind}if _CHECK: return locals()   # p84v2 --check hook" + src[m[0].end():]
        log.append(("--check early-return hook after token construction", 1))
    bad = [p for p in FORBIDDEN if re.search(p, src)]
    assert not bad, f"p82v1 paths or packed TCR/MHC tokens remain after patching: {bad}"
    return src, log


def instantiate_models(ns):
    import torch.nn as nn
    out = []
    for name, obj in list(ns.items()):
        if inspect.isclass(obj) and issubclass(obj, nn.Module) and obj.__module__ in (None, "builtins", "__main__", ns.get("__name__")):
            sig = inspect.signature(obj.__init__)
            kw = {}
            if "mhc_len" in sig.parameters: kw["mhc_len"] = 84
            if "mode" in sig.parameters: kw["mode"] = "struct"
            try:
                mdl = obj(**kw)
            except TypeError as e:
                out.append((name, f"not instantiated without data-specific arguments ({e})")); continue
            shapes = {k: tuple(v.shape) for k, v in list(mdl.named_parameters()) + list(mdl.named_buffers())
                      if any(s in k for s in ("pos", "pair_pos", "sid", "mpos"))}
            out.append((name, f"params {sum(p.numel() for p in mdl.parameters()):,}; {shapes}"))
    return out


def run(consumer, extra=(), cdr3_last_imgt=None, hook_after=None, check_entry=None, argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true", help="train (not used by the rebuild task)")
    a, rest = ap.parse_known_args(argv)
    shim = P.install_grip_shim()
    src0 = Path(consumer).read_text()
    src, log = patch(src0, extra, cdr3_last_imgt, hook_after)
    PATCHED_DIR.mkdir(parents=True, exist_ok=True)
    out = PATCHED_DIR / (Path(consumer).stem + "__p84v2.py")
    out.write_text(src)
    print(f"[{Path(consumer).name} -> p84v2] source sha256 {hashlib.sha256(src0.encode()).hexdigest()[:12]}; "
          f"patched copy {out} ({hashlib.sha256(src.encode()).hexdigest()[:12]})")
    for e, n in log:
        print(f"   edit x{n}: {e}")
    ns = {"__name__": "__main__" if a.run else "p84v2_check", "__file__": consumer, "_P84": P, "_CHECK": not a.run}
    sys.argv = [consumer] + rest
    exec(compile(src, f"{consumer} [p84v2]", "exec"), ns)
    if a.run:
        return None
    res = check_entry(ns) if check_entry else {}
    T, M, pdbs = res["T"], res["M"], [str(p) for p in res["pdbs"]]
    T_ref = P.T_for(pdbs, cdr3_last_imgt=cdr3_last_imgt); M_ref = P.M_for(pdbs)
    print(f"   check: T {T.shape} == canonical: {bool((T == T_ref).all())}; M {M.shape} == canonical: {bool((M == M_ref).all())}; "
          f"complexes {len(pdbs)}; layout SEG_LENS {shim.SEG_LENS}, TCR_LEN {shim.TCR_LEN}")
    for k, v in res.get("arrays", {}).items():
        print(f"   check: {k} shape {tuple(v.shape)}")
    for name, desc in instantiate_models(ns):
        print(f"   model {name}: {desc}")
    ok = bool((T == T_ref).all() and (M == M_ref).all())
    print(f"   RESULT {'OK' if ok else 'MISMATCH'}")
    return ok
