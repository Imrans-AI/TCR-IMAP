#!/usr/bin/env bash
# Progress and ETA for the MD queue -- equilibration and production.
source "$(dirname "$0")/lib.sh"
R0=$MDROOT/md/runs; now=$(date +%s)
fmt(){ local s=$1; printf "%dh %02dm" $((s/3600)) $(((s%3600)/60)); }

eq_done=0; pr_done=0
for pdb in $ORDER; do for rep in 1 2 3; do
  complete "$R0/$pdb/rep$rep/npt"  && eq_done=$((eq_done+1))
  complete "$R0/$pdb/rep$rep/prod" && pr_done=$((pr_done+1))
done; done

echo "=== now running ==="
cur=$(pgrep -fa "gmx mdrun" | grep -o -- "-deffnm [^ ]*" | awk '{print $2}' | head -1)
if [ -n "$cur" ]; then
  st=$(basename "$cur"); rp=$(basename "$(dirname "$cur")"); sy=$(basename "$(dirname "$(dirname "$cur")")")
  t=$(target "$cur"); s=$(livestep "$cur")
  awk -v sy=$sy -v rp=$rp -v st=$st -v s=$s -v t=$t \
      'BEGIN{printf "  %s %s %s : step %d of %d  (%.0f%%)\n", sy, rp, toupper(st), s, t, (t>0?100*s/t:0)}'
else echo "  nothing (GPU idle)"; fi

echo; echo "=== equilibration: $eq_done / 12 replicates complete ==="
if [ -f "$R0/equil.log" ] && [ $eq_done -lt 12 ]; then
  mapfile -t st < <(grep -E "^=== .* equil " "$R0/equil.log" | awk '{print $2}')
  n=${#st[@]}; first=$(date -d "${st[0]}" +%s); last=$(date -d "${st[-1]}" +%s)
  if [ "$n" -ge 2 ]; then avg=$(( (last-first)/(n-1) )); src="measured over $((n-1)) replicates"
  else avg=1500; src="estimate, nothing finished yet"; fi
  cur_el=$(( now-last )); cur_left=$(( avg-cur_el )); [ $cur_left -lt 0 ] && cur_left=0
  left=$(( (12-n)*avg + cur_left ))
  echo "  started          : $(date -d "${st[0]}" '+%H:%M')   (elapsed $(fmt $((now-first))))"
  echo "  per replicate    : $(fmt $avg)   ($src)"
  echo "  remaining        : ~$(fmt $left)"
  echo "  expected finish  : ~$(date -d "@$((now+left))" '+%a %H:%M')"
elif [ $eq_done -eq 12 ]; then echo "  all done -- ready for production (md/prod_start.sh)"; fi

if [ $eq_done -gt 0 ] && compgen -G "$R0/*/rep*/prod.tpr" >/dev/null; then
  echo; echo "=== production: $pr_done / 12 replicates complete ==="
  tot=0
  for pdb in $ORDER; do for rep in 1 2 3; do
    p=$R0/$pdb/rep$rep/prod
    if [ -f "$p.done" ]; then s=50000000
    elif [ "$p" = "$cur" ]; then s=$(livestep "$p")
    else s=$(laststep "$p"); fi
    tot=$((tot+s))
  done; done
  done_ns=$(awk -v s=$tot 'BEGIN{printf "%.1f", s*0.004/1000}')
  rate=146
  if [ -n "$cur" ] && [ "$(basename "$cur")" = prod ]; then
    r=$(awk '/^Started mdrun/{n=0} /Writing checkpoint, step/{a[++n]=$0} END{if(n>=2){print a[n-1]; print a[n]}}' "$cur.log" | awk '
        {step=$4; d=""; for(i=6;i<=NF;i++) d=d" "$i; cmd="date -d \""d"\" +%s"; cmd|getline t; close(cmd)
         if(NR==1){s0=step;t0=t} else if(t>t0) printf "%.1f", (step-s0)*0.004/1000/((t-t0)/86400)}')
    [ -n "$r" ] && rate=$r
  fi
  left_ns=$(awk -v d=$done_ns 'BEGIN{printf "%.1f", 2400-d}')
  echo "  simulated        : $done_ns of 2400 ns"
  echo "  speed            : $rate ns/day"
  echo "  left, continuous : $(awk -v l=$left_ns -v r=$rate 'BEGIN{printf "%.1f days", l/r}')"
fi

echo; bash "$(dirname "$0")/status.sh"
