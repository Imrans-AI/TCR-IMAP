#!/usr/bin/env bash
# Start or resume production in a detached tmux session 'mdprod'. Use each evening.
cd "$(dirname "$0")/.."
pgrep -f run_equil_all.sh >/dev/null && { echo "equilibration still running - not starting production"; exit 1; }
pgrep -f "gmx mdrun" >/dev/null && { echo "an mdrun is already running:"; pgrep -fa "gmx mdrun"; exit 1; }
tmux kill-session -t mdprod 2>/dev/null
tmux new -s mdprod -d
tmux send-keys -t mdprod "source ~/miniconda3/etc/profile.d/conda.sh; conda activate dl_env; cd $PWD; bash md/run_prod_queue.sh 2>&1 | tee -a md/runs/prod_queue.log" C-m
sleep 20; tail -3 md/runs/prod_queue.log; pgrep -fa "gmx mdrun" | cut -c1-120
