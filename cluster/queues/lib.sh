# Queue runner for the gateway: submits the lines of a queue file (resources :: run_* script args) one by one
# whenever a slot is free, records "<jobid> <line>" in the log and retries on the active-job limit.
#   bash cluster/queues/lib.sh <queue file> <log file>
queue=$1 log=$2
cd "$(dirname "$0")/../.."
while IFS= read -r line; do
    [[ -z $line || $line == \#* ]] && continue
    res=${line%% :: *}; args=${line#* :: }
    while :; do
        out=$(slurm submit $res -- $args 2>&1)
        if echo "$out" | grep -q Submitted; then echo "$(echo "$out" | grep -oP 'job \K\d+') $line" >> "$log"; break
        elif echo "$out" | grep -q "active gateway jobs"; then sleep 60
        else echo "FAILED $line :: $out" >> "$log"; break; fi
    done
done < "$queue"
echo "# all submitted $(date +%F_%H:%M)" >> "$log"
