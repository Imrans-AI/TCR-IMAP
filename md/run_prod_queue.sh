#!/usr/bin/env bash
# Production for all 12 replicates, one at a time, interleaved order. Re-run to resume: completed
# replicates are skipped and the current one continues from its checkpoint.
# A crash (mdrun killed by a signal: segfault, abort) is retried from the last checkpoint, at most
# twice per replicate. GROMACS fatal errors (exit 1, e.g. a simulation blow-up) are NOT retried,
# so a genuine instability still stops the queue.
source "$(dirname "$0")/lib.sh"; cd "$MDROOT"
for rep in 1 2 3; do for pdb in $ORDER; do
  tries=0
  while :; do
    echo "=== $(date -Is) prod $pdb rep$rep ==="
    bash md/run_replicate.sh "$pdb" "$rep" prod; rc=$?
    [ $rc -eq 0 ] && break
    if [ $rc -eq 2 ]; then echo "=== $(date -Is) PAUSED $pdb rep$rep (checkpoint saved) ==="; exit 0; fi
    if [ $rc -ge 128 ] && [ $tries -lt 2 ]; then
      tries=$((tries+1))
      echo "=== $(date -Is) CRASH $pdb rep$rep (exit $rc, signal $((rc-128))) - retry $tries/2 from last checkpoint in 60 s ==="
      sleep 60; continue
    fi
    echo "=== $(date -Is) FAILED $pdb rep$rep (exit $rc) ==="; exit $rc
  done
done; done
echo "=== $(date -Is) ALL 12 PRODUCTION RUNS COMPLETE ==="
