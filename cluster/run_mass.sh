#!/bin/bash
# Mass-conservation study: one method of an experiment over the grid sizes of its sweep within [NMIN, NMAX], with
# the irregularity (jitter), interpolation range, weight alpha and MLS order set on the command line (plain
# tokens only, the gateway does not pass quotes/parentheses). CPU, Float64.
#   slurm submit cluster/run_mass.sh cpus=32 mem=32 time=120 -- conservation_burgers 4 epd1 0.4 5.5 1.0 4 1 10000
#   args: NAME ORDER VARIANT JITTER RANGE ALPHA MLS NMIN NMAX
#   VARIANT: plain (no limiter, no MOOD), epd1, epd2 (MOOD U2), vk, mm (minmod, 1-D only); method = RK4MUSCL<ORDER>[(...)]
#            (classical RK4, Rusanov); or the 2026-10-26 setup: ar / ru (SSP stepper, aligned / plain flux, no limiting),
#            ar-epd1, ar-epd2, ar-vk, ar-mm, ru-epd2, ... (see experiments/methods_review.jl; MLS from the MLS argument)
#   MLS: Scheme_MLS_order (0 = coupled to the MOOD order, k = fixed)
set -eo pipefail
source "$(dirname "$0")/common.sh"
name=$1 order=$2 variant=$3 jitter=$4 range=$5 alpha=$6 mls=$7 nmin=$8 nmax=$9
case $variant in
    ar|ru) method=$(review_method "k$order.none.$variant.mls$mls") ;;            # new setup, no limiting
    ar-*|ru-*) method=$(review_method "k$order.${variant#*-}.${variant%%-*}.mls$mls") ;;   # e.g. ar-epd2
    plain) method="RK4MUSCL$order" ;;
    epd1)  method="RK4MUSCL$order(EPD1)" ;;
    epd2)  method="RK4MUSCL$order(EPD2)" ;;
    vk)    method="RK4MUSCL$order(VK)" ;;
    mm)    method="RK4MUSCL$order(minmod)" ;;
    *) echo "unknown variant $variant"; exit 2 ;;
esac
case $name in
    conservation_burgers|Convergence1D|Euler1D)
        dim=1; vary="Grid_Ns__all=filter(n -> $nmin <= n <= $nmax, [ceil(Int, N) for N in 10 .^ (1.5:.1:3.5)])" ;;
    Convergence*)
        dim=2; vary="Grid_Ns__all=filter(n -> $nmin <= n <= $nmax, collect(ceil.(Int, 10 .^ (1.5:.1:3.0))))" ;;
    *) echo "unknown experiment $name"; exit 2 ;;
esac
[ $dim = 1 ] && jit="($jitter,)" || jit="($jitter, $jitter)"
echo "$name / $method  jitter=$jitter R=$range alpha=$alpha MLS=$mls  N in [$nmin, $nmax]  host=$(hostname) threads=$JULIA_NUM_THREADS job=$SLURM_JOB_ID"
export SIM_BACKEND=cpu SIM_T=Float64 SIM_METHODS="$method" SIM_VARY="$vary"
export SIM_SET="Grid_randomness_factor=$jit; Weight_range=$range; Weight_alpha=$alpha; Scheme_MLS_order=$mls"
export SIM_SAVE=${SIM_SAVE:-/work/results} SIM_PROGRESS=${SIM_PROGRESS:-300}
time julia --project=env_headless tools/run_sim.jl "$name"
