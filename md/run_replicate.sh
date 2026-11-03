#!/usr/bin/env bash
# usage: md/run_replicate.sh <pdb> <rep 1|2|3> <equil|prod>
#   equil: minimisation (once per system, shared) -> NVT (seeded) -> NPT, then stop
#   prod : 200 ns production; requires a completed NPT; resumes from its checkpoint
# exit 0 = stage(s) complete | 2 = paused before completion, checkpoint kept | other = error
source "$(dirname "$0")/lib.sh"
set -euo pipefail
PDB=$1; REP=$2; MODE=$3
SEED=${SEEDS[$REP]}
G=$(find "$MDROOT/md/systems/$PDB" -maxdepth 3 -type d -name gromacs | head -1)
P=$MDROOT/md/protocol; EM=$MDROOT/md/runs/$PDB/em; R=$MDROOT/md/runs/$PDB/rep$REP
TOP=$G/topol_hmr.top; NDX=$G/index.ndx; REF=$G/step3_input.gro
GPU="-nb gpu -pme gpu -bonded gpu -update gpu"
export GMX_MAXBACKUP=-1
mkdir -p "$EM" "$R"
free=$(df -BG --output=avail "$MDROOT" | tail -1 | tr -dc 0-9)
[ "$free" -ge 20 ] || { echo "only ${free} GB free - stopping"; exit 3; }
echo "$(date -Is) $MODE pdb=$PDB rep=$REP seed=$SEED git=$(git -C "$MDROOT" rev-parse --short HEAD)" >> "$R/run.info"

stage() {   # stage <prefix> <mdp> <start.gro> <start.cpt or -> <restrained 0|1>
  local pre=$1 mdp=$2 conf=$3 cpt=$4 rst=$5
  complete "$pre" && return 0
  if [ ! -f "$pre.tpr" ]; then
    local args=(-f "$mdp" -c "$conf" -p "$TOP" -n "$NDX" -o "$pre.tpr" -po "${pre}_out.mdp" -maxwarn 0)
    [ "$cpt" != "-" ] && args+=(-t "$cpt")
    [ "$rst" = 1 ] && args+=(-r "$REF")
    gmx grompp "${args[@]}"
  fi
  local rc=0
  gmx mdrun -deffnm "$pre" $GPU -cpi "$pre.cpt" || rc=$?
  complete "$pre" && return 0
  # an interrupted mdrun exits non-zero, so "paused" is read from the log, not the exit code
  if tail -n 300 "$pre.log" | grep -q "Received the .* signal"; then
    echo "$(basename "$pre") paused at step $(laststep "$pre") of $(target "$pre"), checkpoint kept"; exit 2
  fi
  echo "$(basename "$pre") stopped at step $(laststep "$pre") (mdrun exit $rc)"; exit $(( rc ? rc : 1 ))
}

if [ "$MODE" = equil ]; then
  if [ ! -f "$EM/em.done" ]; then
    gmx grompp -f "$G/step4.0_minimization.mdp" -c "$REF" -r "$REF" -p "$TOP" -n "$NDX" \
               -o "$EM/em.tpr" -po "$EM/em_out.mdp" -maxwarn 0
    gmx mdrun -deffnm "$EM/em" -nb gpu
    grep -E "converged" "$EM/em.log" | tail -1 > "$EM/em.done"
  fi
  [ -f "$R/nvt.mdp" ] || sed "s/^gen-seed .*/gen-seed                = $SEED/" "$P/nvt.mdp" > "$R/nvt.mdp"
  stage "$R/nvt" "$R/nvt.mdp" "$EM/em.gro" - 1
  stage "$R/npt" "$P/npt.mdp" "$R/nvt.gro" "$R/nvt.cpt" 1
elif [ "$MODE" = prod ]; then
  complete "$R/npt" || { echo "NPT not complete for $PDB rep$REP - equilibrate first"; exit 4; }
  stage "$R/prod" "$P/prod.mdp" "$R/npt.gro" "$R/npt.cpt" 0
else
  echo "mode must be equil or prod"; exit 5
fi
