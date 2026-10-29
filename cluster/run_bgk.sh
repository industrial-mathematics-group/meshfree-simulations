#!/bin/bash
# BGK on moving grids with particle management (2026-10-28):
#   slurm submit cluster/run_bgk.sh cpus=16 mem=32 time=300 -- CASE METHOD MOVER VOL ADAPT TAU NMIN NMAX
#   CASE   wave (experiments/BGK/bgk_wave.jl, periodic density wave, tools/adapt_study.jl: per-call diagnostics, exact error)
#          sod  (experiments/BGK/bgk_sod.jl, tools/run_sim.jl: L1 against the exact Euler solution, mass drift)
#          sodc (the same Sod through tools/adapt_study.jl: interior conservation defects and events per N)
#          vortex (experiments/BGK/bgk_vortex.jl, 2-D isentropic vortex, one period, tools/adapt_study.jl, N = 50 71 100)
#   optional 9th argument PATCH = on: per-velocity Rusanov dissipation |v_k - w| on moving grids (tools/patches/bgk_ale_flux.jl)
#   METHOD a method label of the config without parentheses-free form: m2 = ARS222MUSCL2, m2vk, m2epd2, m3, m3vk, m3epd2
#   MOVER  none | lagrangian;  VOL shepard | mls;  ADAPT on | off;  TAU relaxation time (e.g. 1e-4, Inf)
set -eo pipefail
source "$(dirname "$0")/common.sh"
case=$1 key=$2 mover=$3 vol=$4 adapt=$5 tau=$6 nmin=$7 nmax=$8 patch=${9:-off}
[ $patch = on ] && export ADAPT_BGK_ALE_FLUX=1
case $key in
    m2) method="ARS222MUSCL2" ;; m2vk) method="ARS222MUSCL2(VK)" ;; m2epd2) method="ARS222MUSCL2(EPD2)" ;;
    m3) method="ARS233MUSCL3" ;; m3vk) method="ARS233MUSCL3(VK)" ;; m3epd2) method="ARS233MUSCL3(EPD2)" ;;
    *) echo "unknown method key $key"; exit 2 ;;
esac
[ $vol = mls ] && v="Grid_volumes=:mls Grid_volumes_degree=2" || v="Grid_volumes=:shepard"
ad=$([ $adapt = on ] && echo true || echo false)
echo "bgk study: case=$case $method mover=$mover vol=$vol adapt=$adapt tau=$tau patch=$patch  host=$(hostname) threads=$JULIA_NUM_THREADS job=$SLURM_JOB_ID"
case $case in
    wave)
        export ADAPT_OUT=/work/results/adapt
        for N in 50 71 100 141 200 283 400 566 800; do
            [ $N -lt $nmin -o $N -gt $nmax ] && continue
            tag="bgkwave_${key}_${mover}_${vol}_${adapt}_tau${tau}_N${N}"
            julia --project=env_headless tools/adapt_study.jl bgk_wave "$method" $tag "Grid_Ns=($N,)" Grid_mover=:$mover PDE_tau=$tau \
                $v Adapt_enabled=$ad || echo "RESULT tag=$tag status=crashed"
        done ;;
    sodc)   # Sod through tools/adapt_study.jl: interior defects (mass, energy; momentum minus 0.9 t) and events per N
        export ADAPT_OUT=/work/results/adapt
        for N in 100 141 200 283 400 566 800 1131 1600; do
            [ $N -lt $nmin -o $N -gt $nmax ] && continue
            tag="bgksod_${key}_${mover}_${vol}_${adapt}_tau${tau}_N${N}"
            julia --project=env_headless tools/adapt_study.jl bgk_sod "$method" $tag "Grid_Ns=($N,)" Grid_mover=:$mover PDE_tau=$tau \
                $v Adapt_enabled=$ad || echo "RESULT tag=$tag status=crashed"
        done ;;
    vortex)
        export ADAPT_OUT=/work/results/adapt
        for N in 50 71 100 141; do
            [ $N -lt $nmin -o $N -gt $nmax ] && continue
            tag="bgkvortex_${key}_${mover}_${vol}_${adapt}_tau${tau}_p${patch}_N${N}"
            julia --project=env_headless tools/adapt_study.jl bgk_vortex "$method" $tag "Grid_Ns=($N,$N)" Grid_mover=:$mover PDE_tau=$tau \
                $v Adapt_enabled=$ad || echo "RESULT tag=$tag status=crashed"
        done ;;
    sod)
        [ $vol = mls ] && vs="Grid_volumes=:mls; Grid_volumes_degree=2" || vs="Grid_volumes=:shepard"
        export SIM_BACKEND=cpu SIM_T=Float64 SIM_METHODS="$method" SIM_SAVE=/work/results SIM_PROGRESS=600
        export SIM_VARY="Grid_Ns__1=filter(n -> $nmin <= n <= $nmax, [100, 141, 200, 283, 400, 566, 800, 1131, 1600]); PDE_tau=[$tau]"
        export SIM_SET="Grid_mover=:$mover; $vs; Adapt_enabled=$ad"
        julia --project=env_headless tools/run_sim.jl bgk_sod ;;
    *) echo "unknown case $case"; exit 2 ;;
esac
