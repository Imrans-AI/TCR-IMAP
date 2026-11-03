# Shared helpers for the MD runners -- sourced, not executed.
command -v gmx >/dev/null 2>&1 || source /usr/local/gromacs-2023.2/bin/GMXRC
declare -A SEEDS=([1]=271828 [2]=314159 [3]=161803)   # fixed; recorded in the criteria document
ORDER="5yxn 7n1f 3vxs 8eo8"                            # alternating strata: H L H L
MDROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

# nsteps the stage was built with, from grompp's processed mdp
target()   { awk -F= '/^ *nsteps/{gsub(/[ \t]/,"",$2); print $2; exit}' "$1_out.mdp" 2>/dev/null || true; }
# last step safely on disk: the log's checkpoint lines, falling back to reading the checkpoint
laststep() {
  local s
  s=$(grep -o "Writing checkpoint, step [0-9]*" "$1.log" 2>/dev/null | tail -1 | awk '{print $4}' || true)
  if [ -z "$s" ] && [ -f "$1.cpt" ]; then
    s=$(gmx dump -cp "$1.cpt" 2>/dev/null | awk '/^ *step *=/{print $3; exit}' || true)
  fi
  echo "${s:-0}"
}
# A stage is complete only if its checkpointed step reached nsteps. An interrupted mdrun still
# writes its final .gro, so the existence of output files is NOT evidence of completion.
complete() {
  [ -f "$1.done" ] && return 0
  local t c; t=$(target "$1"); c=$(laststep "$1")
  if [ -n "$t" ] && [ "$c" -ge "$t" ]; then touch "$1.done"; return 0; fi
  return 1
}

# Current step of a stage, from its energy file. On resume GROMACS truncates the .edr to the
# checkpoint and appends, so its last frame is authoritative; the .log can keep a stale tail
# from a crashed session. Resolution = nstenergy (20 ps in production).
livestep() {
  local t dt
  t=$(gmx check -e "$1.edr" 2>&1 | awk '/Last energy frame read/{print $NF}' || true)
  dt=$(awk -F= '/^ *dt *=/{gsub(/[ \t]/,"",$2); print $2; exit}' "$1_out.mdp" 2>/dev/null || true)
  awk -v t="${t:-0}" -v dt="${dt:-0}" 'BEGIN{ if (dt > 0) printf "%d\n", t/dt + 0.5; else print 0 }'
}
