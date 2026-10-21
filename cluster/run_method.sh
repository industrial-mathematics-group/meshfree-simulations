#!/bin/bash
# One method of an experiment per job (the gateway only passes plain tokens, so methods are short keys):
#   slurm submit cluster/run_method.sh gpus=1 cpus=4 mem=32 time=60 -- Mach3Step m3 600 200 0.4 [cpu|cuda] [Float32|Float64]
# Results go to results/ of the job's run snapshot on the cluster (~/claude-jobs/runs/...).
set -eo pipefail
source "$(dirname "$0")/common.sh"
name=$1 key=$2 nx=$3 ny=$4 cfl=$5 backend=${6:-cuda} T=${7:-Float32}
declare -A METHODS=(
    [m2vk]="RK2MUSCL2(VK)"
    [m2]="RK2MUSCL2(EPD1)"
    [m3]="RK3MUSCL3(EPD1)"            # MLS 0 (coupled to the MOOD order)
    [m3epd2]="RK3MUSCL3(EPD2)"        # MLS 0
    [m3mls3]="RK3MUSCL3(EPD1,MLS3)"
    [m3mls2]="RK3MUSCL3(EPD1,MLS2)"
    [m4]="RK4MUSCL4(EPD1,MLS0)"
    [m4mls4]="RK4MUSCL4(EPD1,MLS4)"
)
method=${METHODS[$key]:?unknown method key '$key' (known: ${!METHODS[*]})}
env=$([ "$backend" = cuda ] && echo env_cuda || echo env_headless)
echo "$name / $method  N=${nx}x${ny} CFL=$cfl  $backend $T  host=$(hostname) threads=$JULIA_NUM_THREADS job=$SLURM_JOB_ID"
export SIM_BACKEND=$backend SIM_T=$T SIM_METHODS="$method" SIM_SET="Grid_Ns=($nx,$ny); Time_CFL=$cfl"
export SIM_SAVE=${SIM_SAVE:-/work/results} SIM_PROGRESS=${SIM_PROGRESS:-120}
time julia --project=$env tools/run_sim.jl "$name"
