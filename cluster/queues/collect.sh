#!/bin/bash
# Waits until every line of a queue has a job id in its log and none of its jobs is queued or running, then collects the
# summary lines of all jobs (RESULT/OK/FAIL/limits, prefixed with the job id) and fetches the *.jld2 files of the
# pipeline jobs into a flat save-path layout per study, <out>/<study>/data/*.jld2 (the layout of tools/*_table.jl):
#   run_method.sh -> euler, run_sweep.sh -> vortex, run_mass.sh -> convergence, run_defect.sh -> defect,
#   run_bgk.sh sod -> bgksod. The tools/adapt_study.jl jobs only report through their RESULT lines.
#     setsid nohup bash cluster/queues/collect.sh <queue> <log> <results.txt> <out dir> &
queue=$1 log=$2 res=$3 out=$4
cd "$(dirname "$0")/../.."
n=$(grep -cvE '^(#|$)' "$queue"); empty=0
while :; do   # done when every queue line has a job id and none of these jobs is still in the gateway queue
    s=$(slurm status 2>&1 < /dev/null)
    if echo "$s" | grep -q "== Queue"; then
        active=$(echo "$s" | sed -n '/== Queue/,/== Recent/p' | grep -oE '^[0-9]+' | grep -cxFf <(grep -oE '^[0-9]+' "$log"))
        [ "$active" -gt 0 ] && empty=0 || empty=$((empty + 1))
    fi
    [ "$(grep -cE '^[0-9]' "$log")" -ge "$n" ] && [ $empty -ge 2 ] && break
    sleep 300
done
echo "# queue drained $(date +%F_%H:%M)" >> "$log"
: > "$res"
grep -E '^[0-9]' "$log" | while read -r j line; do
    slurm log $j 20000 < /dev/null 2>/dev/null |
        grep -E '^RESULT|^OK|^FAIL|TIME LIMIT|CANCELLED|DUE TO|ERROR|bgk study|adapt study|defect study|^\| ' | sed "s/^/$j /" >> "$res"
    set -- ${line#* :: }
    case "$1 $2" in
        cluster/run_method.sh*) study=euler ;; cluster/run_sweep.sh*) study=vortex ;;
        cluster/run_mass.sh*) study=convergence ;; cluster/run_defect.sh*) study=defect ;;
        "cluster/run_bgk.sh sod") study=bgksod ;; *) continue ;;
    esac
    slurm fetch $j < /dev/null > /dev/null 2>&1 || echo "$j fetch failed" >> "$res"
    mkdir -p "$out/$study/data"
    find "data/$j" -name '*.jld2' -exec mv -n {} "$out/$study/data/" \;
    [ -f "data/$j/job-info.txt" ] && mv "data/$j/job-info.txt" "$out/$study/job-info-$j.txt"
    rm -rf "data/$j"
done
echo "# collected $(date +%F_%H:%M)" >> "$log"
