#!/usr/bin/env bash
# Where every replicate stands. Percentages are from the last checkpoint (written every 15 min).
source "$(dirname "$0")/lib.sh"
R0=$MDROOT/md/runs
printf "%-5s %-5s %-6s %-7s %-7s %-7s %s\n" pdb rep EM NVT NPT prod "prod ns"
for pdb in $ORDER; do
  em=$([ -f "$R0/$pdb/em/em.done" ] && echo done || echo -)
  for rep in 1 2 3; do
    cells=""
    for st in nvt npt prod; do
      pre=$R0/$pdb/rep$rep/$st
      if complete "$pre"; then c=done
      elif [ -f "$pre.tpr" ]; then
        # running stage: live step from the log; otherwise the last checkpoint
        if pgrep -f -- "-deffnm $pre " >/dev/null; then st=$(livestep "$pre")
        else st=$(laststep "$pre"); fi
        c=$(awk -v s="$st" -v t="$(target "$pre")" 'BEGIN{printf "%.0f%%",(t>0?100*s/t:0)}')
      else c=-; fi
      cells="$cells$(printf '%-8s' "$c")"
    done
    ns=$(awk -v s="$(laststep "$R0/$pdb/rep$rep/prod")" 'BEGIN{printf "%.1f", s*0.004/1000}')
    printf "%-5s %-5s %-6s %s %s / 200\n" "$pdb" "rep$rep" "$em" "$cells" "$ns"
  done
done
echo; pgrep -fa "gmx mdrun" | cut -c1-110 | sed 's/^/running: /' || echo "running: nothing"
