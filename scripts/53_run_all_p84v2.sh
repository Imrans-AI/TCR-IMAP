#!/usr/bin/env bash
# p84v2 training queue: the eight p84v2 runners with --run, one after another.
# Modelled on md/run_prod_queue.sh. Re-run to resume: a stage whose final outputs all exist is skipped,
# and inside a stage the runners resume from their own results/tables/_partial/p84v2/<stage>/ caches.
# Launch in tmux (the author does this; nothing here was run with --run during setup):
#   tmux new -s p84v2
#   source ~/miniconda3/etc/profile.d/conda.sh; conda activate dl_env; cd ~/PhD_data/TCRMHC_Module_3_TCR-IMAP
#   bash scripts/53_run_all_p84v2.sh 2>&1 | tee -a results/review/train_p84v2_logs/queue_console.log
#
# Flags: --dry-run  print the decisions and exact commands, run nothing, move nothing
#        --force    re-run every stage from scratch: BEFORE any stage runs, every stage's existing p84v2
#                   outputs and caches are MOVED to one backup directory for this invocation,
#                   results/review/train_p84v2_logs/forced_backup_<start time>/<stage>/ (never deleted), and
#                   the three tables derived from them outside the queue are MOVED to .../derived_snapshot/, so
#                   scripts/55_stats_p84v2.py (which writes them only when absent) re-derives them from the new run.
#                   After a failure, resume WITHOUT --force (a second --force would back up the new partial run).
#        --allow-with-md  start even if a gmx mdrun is using the GPU (default: refuse, as md/prod_start.sh does)
#
# SINGLE INSTANCE (added 23 Oct after the queue was launched twice concurrently on 22-23 Oct).
#   flock(1) on results/review/train_p84v2_logs/queue.lock, taken before anything else is checked or written.
#   A second invocation exits at once (exit 4) naming the PID, start time and command of the holder, and
#   touches no output. The kernel drops the lock when the last process holding the descriptor exits,
#   including after kill -9. The runners and pipeline helpers inherit the descriptor ON PURPOSE: if the queue
#   shell is killed while a runner keeps writing, the lock stays held until that runner exits too, so a new
#   queue can never write the same paths as an orphaned runner. queue.lock.info records the holder.
#   Before every stage the queue asserts it still holds the lock and that no queue or runner outside its own
#   process tree is running.
#
# ORDER (justification). No stage reads another stage's output: every runner reads only the p84v2 labels
# (pair_maps_p84v2[_6A].npz), pair_meta_p84v2, tcr_graphs_p84v2 and the frozen metadata/splits. So order is free,
# and is set by the p82v1 wall-clock on this GPU (results/logs/*_p82v1.log), cheapest first, with one exception:
#   1 docking_angle   (16b AngleNet,   99 s)  cheapest: catches an environment/GPU failure in under 2 minutes
#   2 run_all         (24b TCR-IMAP,  ~225 s) the MAIN run (primary tcr70 pairwise result + pairmaps_*_p84v2.npz),
#                     moved ahead of 34b (198 s) so the primary result exists as early as possible and before
#                     everything downstream that depends on it (24c pairwise summaries, Figs 4/6, the MD analysis)
#   3 peptide_aware   (34b,           198 s)
#   4 structcs_setting(31b,           233 s)
#   5 similarity_curve(25b,           293 s)
#   6 teim_res        (36b,           355 s)
#   7 ablations       (29b,           514 s)
#   8 benchmark       (28b,           727 s)  most expensive last (per-residue: Ours, Struct-TCS, Seq-TCS, 5 protocols)
# Expected total ~44 min on the RTX 5070 Ti with the GPU otherwise idle (the 22-23 Oct double run took 78 min).
set -u
ROOT=$(cd "$(dirname "$0")/.." && pwd); cd "$ROOT"
LOGDIR=results/review/train_p84v2_logs; QLOG=$LOGDIR/queue.log
DRY=0; FORCE=0; WITHMD=0
for a in "$@"; do case $a in --dry-run) DRY=1;; --force) FORCE=1;; --allow-with-md) WITHMD=1;;
  *) echo "unknown flag $a"; exit 2;; esac; done

# ---------------- single-instance lock: first thing, before any check or write ----------------
mkdir -p "$LOGDIR"
LOCK=$LOGDIR/queue.lock; LOCKINFO=$LOGDIR/queue.lock.info
exec 9>>"$LOCK"                       # open without truncating; the file itself is never an output
if ! flock -n 9; then
  echo "=== $(date -Is) REFUSED: another p84v2 queue instance holds $LOCK. Nothing was checked, moved or written."
  mine() { local p=$1; while [ -n "$p" ] && [ "$p" -gt 1 ]; do [ "$p" = "$$" ] && return 0
           p=$(ps -o ppid= -p "$p" 2>/dev/null | tr -d ' '); done; return 1; }
  holders=""; for p in $(fuser "$LOCK" 2>/dev/null); do kill -0 "$p" 2>/dev/null && ! mine "$p" && holders="$holders $p"; done
  hp=$(sed -n 's/^pid=\([0-9]*\).*/\1/p' "$LOCKINFO" 2>/dev/null | head -1)
  if [ -n "$hp" ] && case " $holders " in *" $hp "*) true;; *) false;; esac; then
    echo "  running instance: pid $hp, started $(sed -n 's/^started=//p' "$LOCKINFO"), run $(sed -n 's/^run=//p' "$LOCKINFO")"
    echo "  its command: $(sed -n 's/^cmd=//p' "$LOCKINFO")   (cwd $(sed -n 's/^cwd=//p' "$LOCKINFO"))"
  elif [ -n "$hp" ]; then
    echo "  queue.lock.info names pid $hp (started $(sed -n 's/^started=//p' "$LOCKINFO")), which no longer holds the lock:"
    echo "  either that queue was killed and the processes below (its runners) still hold it, or something else took it"
  else
    echo "  no queue.lock.info: the lock was taken by something other than this script"
  fi
  for p in $holders; do
    echo "  lock file open in: $(ps -o pid=,lstart=,args= -p "$p" 2>/dev/null | cut -c1-220)"
  done
  exit 4
fi
STAMP=$(date +%Y%m%dT%H%M%S); RUNID=${STAMP}_$$
CMDLINE=$(tr '\0' ' ' < /proc/$$/cmdline | sed 's/ $//')
STARTED=$(date -Is)
printf 'pid=%s\nstarted=%s\nrun=%s\nmode=%s\ncmd=%s\ncwd=%s\n' "$$" "$STARTED" "$RUNID" \
  "$([ $DRY = 1 ] && echo dry-run || echo run)$([ $FORCE = 1 ] && echo ' --force')" "$CMDLINE" "$ROOT" \
  > "$LOCKINFO.tmp" && mv -f "$LOCKINFO.tmp" "$LOCKINFO"

# is pid $1 this queue, one of its descendants (runners, pipeline helpers), or one of its ancestors (the
# launching shell, e.g. tmux's sh -c '... 53_run_all_p84v2.sh', whose command line names the script)?
ANCESTORS=" "; _p=$PPID
while [ -n "$_p" ] && [ "$_p" -gt 1 ]; do ANCESTORS="$ANCESTORS$_p "; _p=$(ps -o ppid= -p "$_p" 2>/dev/null | tr -d ' '); done
in_my_tree() { local p=$1; case "$ANCESTORS" in *" $p "*) return 0;; esac
               while [ -n "$p" ] && [ "$p" -gt 1 ]; do [ "$p" = "$$" ] && return 0
               p=$(ps -o ppid= -p "$p" 2>/dev/null | tr -d ' '); done; return 1; }
lock_assert() {   # called before every stage: we hold the lock, and nobody else is running the queue or a runner
  grep -q "FLOCK" "/proc/$$/fdinfo/9" 2>/dev/null || { echo "lock lost: fd 9 holds no flock"; return 1; }
  [ "$(sed -n 's/^pid=//p' "$LOCKINFO")" = "$$" ] || { echo "queue.lock.info names another pid"; return 1; }
  local p others=""
  # only processes that EXECUTE the queue or a runner count; one that merely mentions the script name
  # (an editor, less, grep, a wrapper shell's -c string) does not
  for p in $(pgrep -f "^([^ ]*/)?(ba)?sh ([^ ]*/)?scripts/53_run_all_p84v2\.sh|^([^ ]*/)?python[0-9.]* (-u )?([^ ]*/)?scripts/53_[a-z_]+_p84v2\.py --run"); do
    in_my_tree "$p" || others="$others $p"
  done
  others=${others# }
  [ -z "$others" ] || { echo "other queue/runner processes outside this queue:"; ps -o pid=,lstart=,args= -p "${others// /,}" | sed 's/^/      /'; return 1; }
}
qlog() { [ $DRY = 1 ] || echo "$(date -Is) $*" >> "$QLOG"; }

# stage | runner | final outputs (all must exist for the stage to count as complete) | resume cache dirs
STAGES=(
"docking_angle|scripts/53_docking_angle_p84v2.py|results/tables/angle_runs_p84v2.jsonl|results/tables/_partial/p84v2/16b_docking_angle"
"run_all|scripts/53_run_all_p84v2.py|results/tables/protocol_main_p84v2.csv results/tables/pairmaps_random_p84v2.npz results/tables/pairmaps_allele_p84v2.npz results/tables/pairmaps_tcr70_p84v2.npz results/tables/pairmaps_double_p84v2.npz results/tables/pairmaps_temporal_p84v2.npz|results/tables/_partial/p84v2/24b_run_all"
"peptide_aware|scripts/53_peptide_aware_p84v2.py|results/tables/peptide_aware_p84v2.csv results/tables/peptide_preds_p84v2.npz|results/tables/_partial/p84v2/34b_peptide_aware"
"structcs_setting|scripts/53_structcs_setting_p84v2.py|results/tables/structcs_setting_p84v2_4p5A.csv results/tables/structcs_setting_p84v2_6A.csv|results/tables/_partial/p84v2/31b_structcs_45A results/tables/_partial/p84v2/31b_structcs_6A"
"similarity_curve|scripts/53_similarity_curve_p84v2.py|results/tables/similarity_curve_p84v2.csv|results/tables/_partial/p84v2/25b_similarity_curve"
"teim_res|scripts/53_teim_res_p84v2.py|results/tables/teim_res_pairwise_p84v2.npz|"
"ablations|scripts/53_ablations_p84v2.py|results/tables/ablations_p84v2.csv|results/tables/_partial/p84v2/29b_ablations"
"benchmark|scripts/53_benchmark_p84v2.py|results/tables/benchmark_p84v2.csv results/tables/percomplex_random_p84v2.npz results/tables/percomplex_allele_p84v2.npz results/tables/percomplex_tcr70_p84v2.npz results/tables/percomplex_double_p84v2.npz results/tables/percomplex_temporal_p84v2.npz|results/tables/_partial/p84v2/28b_benchmark"
)
# tables computed OUTSIDE the queue from its outputs (scripts/55_stats_p84v2.py, which writes each only when it is
# absent); --force moves them aside with the stage outputs so they are re-derived from the new run, never reused
DERIVED="results/tables/pairwise_ours_p84v2.npz results/tables/percomplex_tcr70_paired_p84v2.tsv results/tables/cluster_macro_p84v2.csv"

say() { echo "=== $(date -Is) $*"; }
die() { say "REFUSED: $*"; exit 3; }

# ---------------- start-up assertions ----------------
fail=0
check() { if eval "$2"; then echo "  ok   $1"; else echo "  FAIL $1"; fail=1; fi; }
say "p84v2 queue start-up checks ($([ $DRY = 1 ] && echo dry-run || echo run)$([ $FORCE = 1 ] && echo ', --force')) pid=$$ run=$RUNID"
echo "  ok   single-instance lock held: $LOCK (pid $$, cmd: $CMDLINE)"
check "conda env dl_env active (CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-unset})" '[ "${CONDA_DEFAULT_ENV:-}" = dl_env ]'
check "GPU visible to PyTorch" 'python -c "import torch,sys; sys.exit(0 if torch.cuda.is_available() else 1)" 2>/dev/null'
FREE_GB=$(df -BG --output=avail "$ROOT" | tail -1 | tr -dc 0-9)
check ">= 20 GB free disk (${FREE_GB} GB)" '[ "$FREE_GB" -ge 20 ]'
hash_ok() { local f=$1; local want; want=$(awk -F'\t' -v f="$f" '$1==f{print $2}' results/review/rebuild_p84v2/artifacts_sha256.tsv)
            [ -f "$f" ] && [ -n "$want" ] && [ "$(sha256sum "$f" | cut -d' ' -f1)" = "$want" ]; }
for f in data/processed/pair_maps_p84v2.npz data/processed/pair_maps_p84v2_6A.npz data/processed/pair_meta_p84v2.parquet data/processed/tcr_graphs_p84v2.npz; do
  check "$f exists and matches the sha256 recorded by the rebuild" "hash_ok $f"
done
if pgrep -f "gmx mdrun" >/dev/null; then
  if [ $WITHMD = 1 ]; then echo "  warn a gmx mdrun is running; starting anyway (--allow-with-md)";
  else echo "  FAIL a gmx mdrun is using the GPU (stop MD with md/prod_stop.sh, or pass --allow-with-md)"; fail=1; fi
else echo "  ok   no gmx mdrun on the GPU"; fi
if lk=$(lock_assert); then echo "  ok   no other queue or runner process outside this one"; else echo "  FAIL $lk"; fail=1; fi
if [ $fail = 1 ]; then [ $DRY = 1 ] && say "dry-run: start-up checks would REFUSE the real run" || die "start-up checks failed"; fi

# ---------------- p82v1 artifacts: snapshot, re-verified after every stage ----------------
SNAP=$(mktemp)
snapshot() { find data/processed results/tables -path '*p82v1*' -type f -printf '%p\t%s\t%T@\n' 2>/dev/null | sort; }
snapshot > "$SNAP"; NP82=$(wc -l < "$SNAP")
echo "  ok   p82v1 artifacts snapshotted: $NP82 files (verified unchanged after each stage)"
RUNLOG=$LOGDIR/run_$RUNID            # this invocation's per-stage logs; earlier logs are never appended to
[ $DRY = 1 ] || mkdir -p "$RUNLOG"
qlog "QUEUE-START run=$RUNID pid=$$ ppid=$PPID user=$(id -un) host=$(hostname) mode=$([ $FORCE = 1 ] && echo force || echo resume) logs=$RUNLOG cmd=[$CMDLINE]"

# ---------------- --force: move every stage's outputs and caches aside BEFORE any stage runs ----------------
if [ $FORCE = 1 ]; then
  BK=$LOGDIR/forced_backup_$STAMP
  say "--force: backup directory for this invocation: $BK"
  nmove=0
  for spec in "${STAGES[@]}"; do
    IFS='|' read -r name runner outputs caches <<< "$spec"
    for o in $outputs $caches; do
      [ -e "$o" ] || continue
      nmove=$((nmove + 1))
      if [ -d "$o" ]; then desc="dir, $(find "$o" -type f | wc -l) files, $(du -sh "$o" | cut -f1)"
      else desc="$(stat -c %s "$o") bytes, sha256 $(sha256sum "$o" | cut -c1-16)"; fi
      if [ $DRY = 1 ]; then echo "  would move  $o  ->  $BK/$name/  ($desc)"
      else mkdir -p "$BK/$name" && mv "$o" "$BK/$name/" || die "could not move $o"
           echo "  moved  $o  ->  $BK/$name/  ($desc)"; qlog "FORCE-MOVE run=$RUNID pid=$$ $o -> $BK/$name/ ($desc)"; fi
    done
  done
  for d in $DERIVED; do
    [ -e "$d" ] || continue
    nmove=$((nmove + 1)); desc="$(stat -c %s "$d") bytes, sha256 $(sha256sum "$d" | cut -c1-16); derived outside the queue, re-derive with 55_stats"
    if [ $DRY = 1 ]; then echo "  would move  $d  ->  $BK/derived_snapshot/  ($desc)"
    else mkdir -p "$BK/derived_snapshot" && mv "$d" "$BK/derived_snapshot/" || die "could not move $d"
         echo "  moved  $d  ->  $BK/derived_snapshot/  ($desc)"; qlog "FORCE-MOVE run=$RUNID pid=$$ $d -> $BK/derived_snapshot/"; fi
  done
  say "--force: $nmove output/cache paths $([ $DRY = 1 ] && echo would be || echo were) moved; every stage will recompute from scratch"
fi

# ---------------- stages ----------------
DONE=(); SKIP=(); FAIL=(); TIMES=()
for spec in "${STAGES[@]}"; do
  IFS='|' read -r name runner outputs caches <<< "$spec"
  for o in $outputs; do case $o in *p84v2*) ;; *) die "stage $name lists a non-p84v2 output $o";; esac; done
  missing=0; for o in $outputs; do [ -e "$o" ] || missing=1; done
  cmd="CUDA_DEVICE_ORDER=PCI_BUS_ID OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 nice -n 15 python -u $runner --run"
  if [ $missing = 0 ] && { [ $FORCE = 0 ] || [ $DRY = 0 ]; }; then
    say "SKIP $name (all outputs exist: $outputs)"; SKIP+=("$name"); qlog "SKIP $name run=$RUNID pid=$$"; continue
  fi
  if [ $DRY = 1 ]; then
    say "WOULD RUN $name: $cmd   > $LOGDIR/run_<run id>/$name.log"; continue
  fi
  lk=$(lock_assert) || { qlog "FAILED $name run=$RUNID pid=$$: lock assertion: $(echo $lk)"; die "lock assertion failed before $name: $lk"; }
  say "START $name ($runner) pid=$$ run=$RUNID"; qlog "START $name run=$RUNID pid=$$ cmd=[$cmd]"
  t0=$(date +%s)
  eval "$cmd" 2>&1 | while IFS= read -r line; do printf '%s %s\n' "$(date -Is)" "$line"; done | tee -a "$RUNLOG/$name.log"
  rc=${PIPESTATUS[0]}; dt=$(( $(date +%s) - t0 ))
  if ! diff -q <(snapshot) "$SNAP" >/dev/null; then
    qlog "FAILED $name run=$RUNID pid=$$: a p82v1 artifact changed"; diff <(snapshot) "$SNAP" | head
    die "p82v1 artifacts changed during stage $name - queue stopped"
  fi
  missing=0; for o in $outputs; do [ -e "$o" ] || missing=1; done
  if [ $rc -ne 0 ] || [ $missing = 1 ]; then
    FAIL+=("$name"); TIMES+=("$name ${dt}s FAILED (exit $rc$([ $missing = 1 ] && echo ', outputs missing'))")
    qlog "FAILED $name run=$RUNID pid=$$ exit=$rc outputs_missing=$missing ${dt}s"
    say "FAILED $name (exit $rc). Queue stopped; later stages not run. Log: $RUNLOG/$name.log"; break
  fi
  for o in $outputs; do qlog "SHA256 $name run=$RUNID $(sha256sum "$o" | cut -d' ' -f1) $(stat -c %s "$o") $o"; done
  DONE+=("$name"); TIMES+=("$name ${dt}s"); qlog "DONE $name ${dt}s run=$RUNID pid=$$"; say "DONE $name (${dt}s)"
done
rm -f "$SNAP"

# ---------------- summary ----------------
say "SUMMARY ($([ $DRY = 1 ] && echo dry-run || echo run)) run=$RUNID pid=$$"
echo "  completed: ${#DONE[@]}  ${DONE[*]:-}"
echo "  skipped:   ${#SKIP[@]}  ${SKIP[*]:-}"
echo "  failed:    ${#FAIL[@]}  ${FAIL[*]:-}"
for t in "${TIMES[@]:-}"; do [ -n "$t" ] && echo "  wall-clock $t"; done
qlog "QUEUE-END run=$RUNID pid=$$ completed=${#DONE[@]} skipped=${#SKIP[@]} failed=${#FAIL[@]}"
[ ${#FAIL[@]} -eq 0 ] || exit 1
