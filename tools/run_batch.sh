#!/usr/bin/env bash
# Runs every method of an experiment in its own Julia process (device runs in one process can corrupt each
# other's buffers with AMDGPU 2.9.1, see HyperCloudCore docs/src/gpu_roadmap.md) and reruns a failed method
# once with SIM_FORCE=1 to tell a reproducible (numerical) failure from a flaky one.
#
#     tools/run_batch.sh NAME "method1;method2;..."      (other SIM_* variables are passed through)
set -uo pipefail
cd "$(dirname "$0")/.."
name=$1; IFS=';' read -ra methods <<< "$2"
log="logs/${name}_$(date +%Y%m%d-%H%M%S).log"
mkdir -p logs
for m in "${methods[@]}"; do
    echo "=== $name / $m  (T=${SIM_T:-Float32}, ${SIM_BACKEND:-gpu}, SET='${SIM_SET:-}')" | tee -a "$log"
    SIM_METHODS="$m" julia -t "${SIM_THREADS:-6}" --project="${SIM_ENV:-env_gpu}" tools/run_sim.jl "$name" >> "$log" 2>&1
    rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "--- failed (rc=$rc), rerunning once" | tee -a "$log"
        SIM_FORCE=1 SIM_METHODS="$m" julia -t "${SIM_THREADS:-6}" --project="${SIM_ENV:-env_gpu}" tools/run_sim.jl "$name" >> "$log" 2>&1
        rc=$?
        [[ $rc -ne 0 ]] && echo "--- failed again (rc=$rc): reproducible" | tee -a "$log" || echo "--- second attempt OK: flaky first run" | tee -a "$log"
    fi
    grep -E "^(OK|FAIL) " "$log" | tail -n 1
done
echo "log: $log"
