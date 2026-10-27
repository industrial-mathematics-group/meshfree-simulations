#!/bin/bash
# Particle-management study (tools/adapt_study.jl): one case / method / setting over a range of grid sizes per job.
#   slurm submit cluster/run_adapt.sh cpus=16 mem=32 time=300 -- CASE K VOL ADAPT MAXGAP MINGAP SPD EXCL NMIN NMAX
#   CASE   b1 (1-D Burgers, periodic Gaussian, Lagrangian mover, t = 2, after the shock)
#          b1pre (same, t = 0.5, before the shock), lin (1-D advection, prescribed sinusoidal mover A = 0.5, t = 2)
#          b2 (2-D Burgers, periodic Gaussian, Lagrangian mover, t = 1.5, after the shock)
#          vortex (isentropic vortex, Lagrangian mover, one period t = 10)
#   K      scheme order: RK2MUSCL2 / RK3MUSCL3 / RK4SSPMUSCL4 with EPD2, AlignedRusanov, frozen MLS 2; 0 = MUSCL2 without MOOD
#   VOL    shepard | mls (volume estimator; mls = degree 2)
#   ADAPT  on | off;  MAXGAP MINGAP SPD (samples_per_dx) EXCL (exclude_reduced_order: 1|0) the adaptivity settings
#   NMIN NMAX: the grid sizes of the case's list within this range
# Writes results/adapt/<tag>.{csv,jld2}; fetch with `slurm fetch <id> --pattern '*.csv'` for the step logs.
set -eo pipefail
source "$(dirname "$0")/common.sh"
case=$1 K=$2 vol=$3 adapt=$4 maxgap=$5 mingap=$6 spd=$7 excl=$8 nmin=$9 nmax=${10}
case $K in 0) method="RK2MUSCL2(none,AR,MLS2)" ;; *) method=$(review_method "k$K.epd2.ar.mls2") ;; esac
common="Grid_randomness_factor=(0.2,) Grid_mover=:lagrangian"
case $case in
    b1)     name=Convergence1D;  Ns="100 141 200 283 400 566 800 1131 1600"; extra="PDE_name=:burgers tmax=2.0" ;;
    b1pre)  name=Convergence1D;  Ns="100 141 200 283 400 566 800 1131 1600"; extra="PDE_name=:burgers tmax=0.5" ;;
    lin)    name=Convergence1D;  Ns="100 141 200 283 400 566 800 1131 1600"; extra="tmax=2.0"
            common="Grid_randomness_factor=(0.2,) Grid_mover=:sinusoidal Grid_mover_amplitude=0.5" ;;
    b2)     name=ConvergenceBurgers; Ns="50 71 100 141 200 283"; extra="tmax=1.5"
            common="Grid_randomness_factor=(0.2,0.2) Grid_mover=:lagrangian" ;;
    vortex) name=IsentropicVortex; Ns="50 71 100 141 200"; extra="tmax=10.0"
            common="Grid_randomness_factor=(0.2,0.2) Grid_mover=:lagrangian" ;;
    *) echo "unknown case $case"; exit 2 ;;
esac
[ $vol = mls ] && v="Grid_volumes=:mls Grid_volumes_degree=2" || v="Grid_volumes=:shepard"
a="Adapt_enabled=$([ $adapt = on ] && echo true || echo false) Adapt_max_gap=$maxgap Adapt_min_gap=$mingap Adapt_samples_per_dx=$spd Adapt_exclude_reduced_order=$([ $excl = 1 ] && echo true || echo false)"
export ADAPT_OUT=/work/results/adapt
echo "adapt study: case=$case K=$K ($method) vol=$vol adapt=$adapt gaps=$maxgap/$mingap spd=$spd excl=$excl  host=$(hostname) threads=$JULIA_NUM_THREADS job=$SLURM_JOB_ID"
for N in $Ns; do
    [ $N -lt $nmin -o $N -gt $nmax ] && continue
    case $name in Convergence1D) g="Grid_Ns=($N,)" ;; *) g="Grid_Ns=($N,$N)" ;; esac
    tag="${case}_K${K}_${vol}_${adapt}_g${maxgap}-${mingap}_s${spd}_e${excl}_N${N}"
    julia --project=env_headless tools/adapt_study.jl $name "$method" $tag $g $common $extra $v $a || echo "RESULT tag=$tag status=crashed"
done
