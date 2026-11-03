#!/usr/bin/env bash
# Pause production cleanly. Sends ONE interrupt to mdrun only: it stops at the next neighbour-
# search step, writes a checkpoint and exits. Never press Ctrl-C in the tmux window instead.
pid=$(pgrep -f "mdrun -deffnm .*/md/runs/.*/prod" | head -1)
[ -n "$pid" ] || { tmux kill-session -t mdprod 2>/dev/null && echo "no mdrun running - stopped the idle queue" || echo "no production running"; exit 0; }
echo "interrupting mdrun (pid $pid), waiting for checkpoint..."; kill -INT "$pid"
for i in $(seq 1 120); do kill -0 "$pid" 2>/dev/null || break; sleep 2; done
kill -0 "$pid" 2>/dev/null && { echo "still running after 4 min - check manually"; exit 1; }
sleep 3; tmux kill-session -t mdprod 2>/dev/null
echo "paused cleanly."; bash "$(dirname "$0")/status.sh"
