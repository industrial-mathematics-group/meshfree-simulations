#!/bin/bash
# 2-D BGK runs of study/bgk-lagrangian-contact (experiments/BGK/bgk_sod2d.jl), one method and resolution per job:
#   slurm submit cluster/run_bgk2d.sh gpus=1 cpus=16 mem=32 time=300 -- CONFIG KEY NX NY [patch]
#   KEY  <vk|epd2|ars233>-<s|l>: ARS222MUSCL2(VK), ARS222MUSCL2(EPD2), ARS233MUSCL3(EPD2) on the static (s) or
#        Lagrangian (l) grid (method names "<method> static" / "<method> lagrangian" of the config)
#   patch  loads tools/patches/bgk_ale_flux.jl (per-velocity dissipation |v_k − w| on moving grids); the results then go
#          to results/alepatch/ (same parameters, so they must not mix with the unpatched runs)
# CPU backend, Float64 (the cores of whatever node the job gets). Results: results/*.jld2 of the run snapshot.
set -eo pipefail
source "$(dirname "$0")/common.sh"
config=$1 key=$2 nx=$3 ny=$4 patch=${5:-off}
case ${key%-*} in
    vk) m="ARS222MUSCL2(VK)" ;; epd2) m="ARS222MUSCL2(EPD2)" ;; ars233) m="ARS233MUSCL3(EPD2)" ;;
    *) echo "unknown method key $key"; exit 2 ;;
esac
case ${key#*-} in s) m="$m static" ;; l) m="$m lagrangian" ;; *) echo "unknown grid in $key"; exit 2 ;; esac
export SIM_BACKEND=cpu SIM_T=Float64 SIM_METHODS="$m" SIM_SET="Grid_Ns=($nx,$ny)" SIM_PROGRESS=600 SIM_SAVE=/work/results
if [ "$patch" = patch ]; then export SIM_PATCHES=tools/patches/bgk_ale_flux.jl SIM_SAVE=/work/results/alepatch; fi
echo "bgk2d: $config / $m N=${nx}x${ny} patch=$patch host=$(hostname) threads=$JULIA_NUM_THREADS job=$SLURM_JOB_ID"
julia --project=env_headless -e 'using Pkg; Pkg.status(["HyperCloud", "HyperCloudCore", "PDEStudioCore"])'
time julia --project=env_headless tools/run_sim.jl "$config"
