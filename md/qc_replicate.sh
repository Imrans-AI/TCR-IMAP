#!/usr/bin/env bash
# Integrity QC of one finished production replicate. Read-only on the run; scratch work in /tmp.
# Deliberately does NOT compute per-pair occupancy or anything involving TCR-IMAP scores -- that is
# the predeclared outcome analysis (criteria section 6), run only once all 12 replicates exist.
# usage: md/qc_replicate.sh <pdb> <rep>
source "$(dirname "$0")/lib.sh"
PDB=$1; REP=$2
P=$MDROOT/md/runs/$PDB/rep$REP/prod
G=$(find "$MDROOT/md/systems/$PDB" -maxdepth 3 -type d -name gromacs | head -1)
T=$(mktemp -d /tmp/qc_${PDB}_rep${REP}_XXXX); cd "$T"; export GMX_MAXBACKUP=-1
echo "################ $PDB rep$REP ################"

cpt=$(gmx dump -cp "$P.cpt" 2>/dev/null | awk '/^ *step *=/{print $3; exit}')
echo "  checkpoint step : $cpt / 50000000  $([ "$cpt" = 50000000 ] && echo OK || echo '<< INCOMPLETE')"
ev=$(grep -oE "(CRASH|FAILED|PAUSED) $PDB rep$REP" "$MDROOT/md/runs/prod_queue.log" | awk '{print $1}' | sort | uniq -c | xargs)
echo "  queue events    : ${ev:-none}   | mdrun sessions $(grep -c '^Started mdrun' "$P.log") | LINCS warnings $(grep -c 'LINCS WARNING' "$P.log")"

ck=$(gmx check -f "$P.xtc" 2>&1)
fr=$(echo "$ck" | awk '/^Coords/{print $2}'); dt=$(echo "$ck" | awk '/^Coords/{print $3}')
last=$(echo "$ck" | tr '\r' '\n' | awk '/Last frame/{print $NF}'); irr=$(echo "$ck" | grep -c "don't match")
echo "  trajectory      : $fr frames every $dt ps, last t=$last ps, irregular steps $irr  $([ "$fr" = 10001 ] && [ "$irr" = 0 ] && echo OK || echo '<< CHECK')"
echo "  energy file     : last t=$(gmx check -e "$P.edr" 2>&1 | awk '/Last energy frame read/{print $NF}') ps"

printf "Temperature\nPressure\nDensity\n" | gmx energy -f "$P.edr" -b 20000 2>/dev/null \
  | awk '/^(Temperature|Pressure|Density)/{printf "  %-15s : avg %9.2f  sd %7.2f  drift(20-200 ns) %+8.2f\n",$1,$2,$4,$5}'

echo SOLU | gmx convert-tpr -s "$P.tpr" -n "$G/index.ndx" -o prot.tpr >/dev/null 2>&1
echo 0 | gmx trjconv -s prot.tpr -f "$P.xtc" -o w.xtc  -pbc whole  -dt 1000 >/dev/null 2>&1
echo 0 | gmx trjconv -s prot.tpr -f w.xtc     -o nj.xtc -pbc nojump          >/dev/null 2>&1
echo 0 | gmx trjconv -s prot.tpr -f nj.xtc    -o frames.gro -dt 10000        >/dev/null 2>&1

python - frames.gro "$G" <<'PY'
import sys, re, numpy as np
from scipy.spatial import cKDTree
gro, G = sys.argv[1], sys.argv[2]
def frames(p):
    L = open(p).read().splitlines(); i = 0
    while i < len(L) and L[i].strip():
        n = int(L[i+1]); a = L[i+2:i+2+n]
        m = re.search(r"t=\s*([0-9.eE+-]+)", L[i]); t = float(m.group(1)) if m else 0.0
        yield t, [x[10:15].strip() for x in a], np.array([[float(x[20:28]), float(x[28:36]), float(x[36:44])] for x in a])
        i += n + 3
sizes = []
blk = re.search(r"\[\s*molecules\s*\](.*)", open(f"{G}/topol.top").read(), re.S).group(1)
for l in blk.splitlines():
    t = l.split(";")[0].split()
    if len(t) < 2 or not t[0].startswith("PRO"): continue
    on, c = False, 0
    for x in open(f"{G}/toppar/{t[0]}.itp"):
        s = x.strip()
        if s.startswith("["): on = s.replace(" ", "") == "[atoms]"; continue
        if on and s and not s.startswith(";"): c += 1
    sizes.append((t[0], c))
role = dict(zip([n for n, _ in sorted(sizes, key=lambda x: x[1])], ["pep", "b2m", "TCRa", "TCRb", "MHC"]))
idx, st = {}, 0
for n, c in sizes: idx[role[n]] = np.arange(st, st + c); st += c
PAIRS = [("pep","MHC"), ("pep","TCRa"), ("pep","TCRb"), ("TCRa","MHC"), ("TCRb","MHC"), ("MHC","b2m"), ("TCRa","TCRb")]
def kabsch(A, B):
    A = A - A.mean(0); B = B - B.mean(0); U, S, Vt = np.linalg.svd(A.T @ B)
    D = np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))]); return np.sqrt((((A @ (U @ D @ Vt)) - B) ** 2).sum(1).mean())
ts, rms, cont, ref = [], [], {p: [] for p in PAIRS}, None
for t, names, x in frames(gro):
    assert len(x) == st, f"atom count {len(x)} != {st}"
    heavy = np.array([not nm.startswith("H") for nm in names]); bb = np.array([nm in ("N","CA","C") for nm in names])
    if ref is None: ref = x[bb]
    ts.append(t); rms.append(kabsch(ref, x[bb]))
    for a, b in PAIRS:
        ia, ib = idx[a][heavy[idx[a]]], idx[b][heavy[idx[b]]]
        cont[(a, b)].append(int(cKDTree(x[ia]).count_neighbors(cKDTree(x[ib]), 0.45)))
ts, rms = np.array(ts) / 1000, np.array(rms)
late = rms[ts >= 100]
print(f"  backbone RMSD   : mean(100-200 ns) {late.mean():.2f} nm, max {rms.max():.2f}, final {rms[-1]:.2f}  ({len(ts)} frames, 10 ns apart)")
print("  interfaces      : heavy-atom contacts < 4.5 A, start -> min -> end")
for p in PAIRS:
    c = cont[p]; flag = "  << CHECK" if c[0] > 20 and (c[-1] < 0.25 * c[0] or min(c) < 0.10 * c[0]) else ""
    print(f"     {p[0]+'-'+p[1]:10s} {c[0]:5d} -> {min(c):5d} -> {c[-1]:5d}{flag}")
PY
cd / && rm -rf "$T"
