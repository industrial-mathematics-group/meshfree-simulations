#!/bin/bash
# Conservation-defect experiments (experiments/ConservationDefect) through the PDEStudioCore pipeline, one or more
# configs per job; the tables go to the log (and results/defect/tables_<configs>.md), the runs to results/defect/.
#   slurm submit cluster/run_defect.sh cpus=32 mem=64 time=120 -- defect_shock [defect_floor ...]
set -eo pipefail
source "$(dirname "$0")/common.sh"
echo "defect study: $*  host=$(hostname) threads=$JULIA_NUM_THREADS job=$SLURM_JOB_ID"
mkdir -p results/defect
export DEFECT_SAVE=/work/results/defect DEFECT_TABLES="results/defect/tables_$(IFS=_; echo "$*").md"
julia --project=env_headless tools/run_defect.jl "$@"
