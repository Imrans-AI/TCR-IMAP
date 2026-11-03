#!/usr/bin/env bash
# Equilibrate all 12 replicates: minimisation once per system, then NVT + NPT per replicate (~5 h).
source "$(dirname "$0")/lib.sh"; cd "$MDROOT"
for rep in 1 2 3; do for pdb in $ORDER; do
  echo "=== $(date -Is) equil $pdb rep$rep ==="
  bash md/run_replicate.sh "$pdb" "$rep" equil || { rc=$?; echo "STOPPED $pdb rep$rep (exit $rc)"; exit $rc; }
done; done
echo "=== $(date -Is) ALL 12 EQUILIBRATIONS COMPLETE ==="
