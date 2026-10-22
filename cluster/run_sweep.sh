#!/bin/bash
# Part of a convergence sweep per job: one method, the grid sizes of the experiment's sweep between NMIN and NMAX
# (the gateway only passes plain tokens, so the sweep is restricted by a range instead of a list):
#   slurm submit cluster/run_sweep.sh cpus=32 mem=64 time=300 -- Convergence1D RK4MUSCL4 1 10000 cpu Float64
#   slurm submit cluster/run_sweep.sh gpus=1 gputype=v100 cpus=4 mem=64 time=300 -- ConvergenceAdvection RK4MUSCL5 300 1000 cuda Float64
#   ... -- ConvergenceAdvection RK4MUSCL4 251 1000 cuda Float64 jitter0     (regular lattice instead of jitter 0.2)
# The N values are those of the experiment files (experiments/Convergence/*.jl), filtered to [NMIN, NMAX]; the
# tasks are the same as a full `@sim` of the experiment, so the results can be combined.
set -eo pipefail
source "$(dirname "$0")/common.sh"
name=$1 method=$2 nmin=$3 nmax=$4 backend=${5:-cpu} T=${6:-Float64} opt=$7   # opt: jitter0 = regular lattice
case $name in
    Convergence1D) vary="Grid_Ns__1=filter(n -> $nmin <= n <= $nmax, [ceil(Int, N) for N in 10 .^ (1.5:.1:3.5)])" ;;
    Convergence*)  vary="Grid_Ns__all=filter(n -> $nmin <= n <= $nmax, collect(ceil.(Int, 10 .^ (1.5:.1:3.0))))" ;;
    *) echo "unknown sweep experiment $name"; exit 2 ;;
esac
env=$([ "$backend" = cuda ] && echo env_cuda || echo env_headless)
echo "$name / $method  N in [$nmin, $nmax]  $backend $T  host=$(hostname) threads=$JULIA_NUM_THREADS job=$SLURM_JOB_ID"
export SIM_BACKEND=$backend SIM_T=$T SIM_METHODS="$method" SIM_VARY="$vary"
[ "$opt" = jitter0 ] && export SIM_SET="Grid_randomness_factor=(0.0, 0.0)"
export SIM_SAVE=${SIM_SAVE:-/work/results} SIM_PROGRESS=${SIM_PROGRESS:-300}
time julia --project=$env tools/run_sim.jl "$name"
