#!/usr/bin/env bash
cd /home/node/.claude/jobs/a4b4be7c/tmp
for cfl in 0.4 0.8; do for K in 3 4; do for s in 0 1; do
    echo -n "CFL=$cfl "; CFL=$cfl julia -t 8 --project=env_ssp ssp_study.jl vac $K 400 $s 2>&1 | tail -1
done; done; done
echo VAC DONE
