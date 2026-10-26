#!/bin/bash
# Submits the 2026-10-26 queues in order (convergence, mass, Euler, vortex), then waits until the gateway queue is
# empty. Logs: cluster/queues/logs/<queue>.log (job id + line). Run with nohup/background; safe to re-run: it skips
# queue lines that already have a job id in the log.
cd "$(dirname "$0")/../.."; mkdir -p cluster/queues/logs
for q in convergence mass euler vortex; do
    queue=cluster/queues/2026-10-26-$q.txt; log=cluster/queues/logs/2026-10-26-$q.log; touch "$log"
    tmp=$(mktemp); grep -vxFf <(sed -E 's/^[0-9]+ //' "$log") "$queue" > "$tmp" || true
    bash cluster/queues/lib.sh "$tmp" "$log"; rm -f "$tmp"
done
while :; do s=$(slurm status 2>&1); echo "$s" | grep -q "== Queue" || { sleep 60; continue; }
    echo "$s" | sed -n '/== Queue/,/== Recent/p' | grep -qE '^[0-9]' || break; sleep 120; done
echo "# queue drained $(date +%F_%H:%M)" >> cluster/queues/logs/2026-10-26-vortex.log
