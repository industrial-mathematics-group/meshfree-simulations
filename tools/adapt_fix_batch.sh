#!/usr/bin/env bash
# Local before/after rerun of the particle-management study (tools/adapt_study.jl, cases as cluster/run_adapt.sh)
# against two HyperCloudCore checkouts, e.g. main and claude/adaptivity-conservation:
#
#   ENV_MAIN=/path/to/env_main ENV_FIX=/path/to/env_fix bash tools/adapt_fix_batch.sh
#
# Each environment has the registered HyperCloud plus `Pkg.develop(path = <HyperCloudCore checkout>)` (and
# StaticArrays, Logging); create them outside the repo. RESULT lines (prefixed with env=env_main|env_fix) go to
# results/adapt_fix/all.log; tables: python3 tools/adapt_fix_table.py results/adapt_fix/all.log
# Every run compiles anew (about 1–2 min each); the whole batch took about 70 min on 8 threads.
set -o pipefail
cd "$(dirname "$0")/.."
OUT=results/adapt_fix; LOG=$OUT/all.log
mkdir -p $OUT; : > $LOG

run() {   # run ENVNAME CASE K VOL ADAPT N...
    local envname=$1 case=$2 K=$3 vol=$4 adapt=$5; shift 5
    local envdir; [ $envname = env_main ] && envdir=$ENV_MAIN || envdir=$ENV_FIX
    local method name extra common v a g tag
    case $K in 2) method="RK2MUSCL2(EPD2,AR,MLS2)" ;; 4) method="RK4SSPMUSCL4(EPD2,AR,MLS2)" ;; esac
    common="Grid_randomness_factor=(0.2,) Grid_mover=:lagrangian"
    case $case in
        b1)     name=Convergence1D;      extra="PDE_name=:burgers tmax=2.0" ;;
        b1pre)  name=Convergence1D;      extra="PDE_name=:burgers tmax=0.5" ;;
        lin)    name=Convergence1D;      extra="tmax=2.0"; common="Grid_randomness_factor=(0.2,) Grid_mover=:sinusoidal Grid_mover_amplitude=0.5" ;;
        b2)     name=ConvergenceBurgers; extra="tmax=1.5"; common="Grid_randomness_factor=(0.2,0.2) Grid_mover=:lagrangian" ;;
        vortex) name=IsentropicVortex;   extra="tmax=10.0"; common="Grid_randomness_factor=(0.2,0.2) Grid_mover=:lagrangian" ;;
    esac
    [ "$vol" = mls ] && v="Grid_volumes=:mls Grid_volumes_degree=2" || v="Grid_volumes=:shepard"
    a="Adapt_enabled=$([ "$adapt" = on ] && echo true || echo false) Adapt_max_gap=1.6 Adapt_min_gap=0.5 Adapt_samples_per_dx=2 Adapt_exclude_reduced_order=true"
    for N in "$@"; do
        case $name in Convergence1D) g="Grid_Ns=($N,)" ;; *) g="Grid_Ns=($N,$N)" ;; esac
        tag="${case}_K${K}_${vol}_${adapt}_N${N}"
        echo "== $envname $tag" >> $LOG
        ADAPT_OUT=$OUT/$envname timeout 30m julia -t 8 --project=$envdir tools/adapt_study.jl $name "$method" $tag $g $common $extra $v $a 2>&1 \
            | grep "^RESULT" | sed "s/^RESULT /RESULT env=$envname /" >> $LOG || echo "RESULT env=$envname tag=$tag status=crashed" >> $LOG
    done
}

for e in env_main env_fix; do
    run $e lin    2 mls on  200 400 800 1600
    run $e b1pre  2 mls on  200 400 800 1600
    run $e b1     2 mls on  200 400 800
    run $e vortex 2 shepard on 71
    run $e vortex 2 mls on  50 71
    run $e b2     2 mls on  71 100
done
# Adaptivity off: untouched by the library changes (the driver sets the stepper's estimator itself)
run env_fix lin    2 mls off 200 400 800 1600
run env_fix b1pre  2 mls off 200 400 800 1600
run env_fix b1     2 mls off 200 400 800
run env_fix vortex 2 mls off 50 71
run env_fix b2     2 mls off 71 100
echo "BATCH DONE" >> $LOG
