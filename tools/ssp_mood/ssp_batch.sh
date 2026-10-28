#!/usr/bin/env bash
# SSP substep checks (on/off) for K = 2, 3, 4 and the K = 1 orientation (fix vs main).
cd /home/node/.claude/jobs/a4b4be7c/tmp
R="julia -t 8 --project=env_ssp ssp_study.jl"
for K in 2 3 4; do
    for s in 0 1; do
        $R sod  $K 400 $s 2>&1 | tail -1
        $R r123 $K 400 $s 2>&1 | tail -1
        $R box  $K 100 $s 2>&1 | tail -1
    done
done
for K in 3 4; do
    for s in 0 1; do $R vortex $K 141 $s 1e-2 2>&1 | tail -1; done
done
for c in "sod 400" "r123 400" "box 100"; do
    set -- $c
    echo -n "fix  "; julia -t 8 --project=env_ssp ssp_study.jl $1 1 $2 1 2>&1 | tail -1
    echo -n "main "; julia -t 8 --project=env_main2 ssp_study.jl $1 1 $2 1 2>&1 | tail -1
done
echo SSP DONE
