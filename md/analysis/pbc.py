"""Periodic-boundary handling for the characterisation measures (criteria document, "Characterisation recomputed with
periodic-boundary handling", 24 October 2026).

`whole_nojump_xtc` makes a production trajectory whole and jump-free by exactly the procedure of md/qc_replicate.sh:
convert the run input to the protein-only SOLU group, then `gmx trjconv -pbc whole`, then `gmx trjconv -pbc nojump`.
The only difference is that every frame is kept (qc_replicate.sh subsamples with -dt 1000).

Only the coordinate-based characterisation measures need this. Occupancy, the docked-frame counts and the hydrogen
bonds use periodic-aware distances and read the original prod.xtc.
"""
import subprocess
from pathlib import Path

GMXRC = "/usr/local/gromacs-2023.2/bin/GMXRC"


def _sh(cmd, cwd):
    full = f"command -v gmx >/dev/null 2>&1 || source {GMXRC}; export GMX_MAXBACKUP=-1; {cmd}"
    r = subprocess.run(["bash", "-c", full], cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"failed ({r.returncode}): {cmd}\n{r.stderr[-2000:]}")
    return r.stdout + r.stderr


def xtc_frames(xtc, cwd):
    """(frame count, first time, last time) of an xtc, from `gmx check`."""
    out = _sh(f"gmx check -f '{xtc}' 2>&1", cwd)
    n = [l.split()[1] for l in out.splitlines() if l.startswith("Coords")]
    first = [l for l in out.splitlines() if "Reading frame       0 time" in l]
    last = [x for x in out.replace("\r", "\n").splitlines() if "Last frame" in x]
    return int(n[0]) if n else -1, first[0].split()[-1] if first else None, last[0].split()[-1] if last else None


def whole_nojump_xtc(pdb, rep, out_dir, root):
    """Write <out_dir>/<pdb>_rep<rep>/nj.xtc: prod.xtc made whole, then jump-free, every frame. Returns its path.
    The intermediate whole-only trajectory is deleted. Atom order and count are those of the protein-only SOLU group,
    i.e. of prod.xtc itself, so every existing atom selection applies unchanged."""
    import glob
    P = Path(root) / f"md/runs/{pdb}/rep{rep}/prod"
    G = Path(glob.glob(str(Path(root) / f"md/systems/{pdb}/charmm-gui-*/gromacs"))[0])
    d = Path(out_dir) / f"{pdb}_rep{rep}"
    d.mkdir(parents=True, exist_ok=True)
    _sh(f"echo SOLU | gmx convert-tpr -s '{P}.tpr' -n '{G}/index.ndx' -o prot.tpr", d)
    _sh(f"echo 0 | gmx trjconv -s prot.tpr -f '{P}.xtc' -o w.xtc -pbc whole", d)
    _sh("echo 0 | gmx trjconv -s prot.tpr -f w.xtc -o nj.xtc -pbc nojump", d)
    (d / "w.xtc").unlink()
    return d / "nj.xtc"
