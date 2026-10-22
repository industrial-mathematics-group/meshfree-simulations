#!/bin/bash
# Smoke test of a backend on the cluster: Euler2D on a 100x100 grid with two methods.
#   slurm submit cluster/backend_test.sh cpus=8 -- cpu
#   slurm submit cluster/backend_test.sh gpus=1 cpus=4 -- cuda
set -eo pipefail
source "$(dirname "$0")/common.sh"
backend=${1:-cpu}
env=$([ "$backend" = cuda ] && echo env_cuda || echo env_headless)
echo "backend=$backend env=$env host=$(hostname) threads=$JULIA_NUM_THREADS"
if [ "$backend" = cuda ]; then
    nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
    julia --project=$env -e 'using CUDACore; println("CUDA runtime ", CUDACore.runtime_version(), ", driver ", CUDACore.driver_version(), ", device ", CUDACore.name(CUDACore.device()), " (", CUDACore.capability(CUDACore.device()), ")")'
fi
export SIM_BACKEND=$backend SIM_SET="Grid_Ns=(100,100)" SIM_METHODS="RK2MUSCL2(EPD1);RK3MUSCL3(EPD1)"
export SIM_SAVE=/tmp/sim-results-$SLURM_JOB_ID SIM_PROGRESS=10
time julia --project=$env tools/run_sim.jl Euler2D
