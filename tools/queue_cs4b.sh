#!/usr/bin/env bash
# clainshock4 with the PAD + divergence-order fix (worktrees, env_pad), GPU Float32, configured 600x600 grid
cd "$(dirname "$0")/.."
export SIM_ENV=env_pad SIM_SAVE=results_pad SIM_PROGRESS=300
for cfl in 0.2 0.4; do
  SIM_SET="Time_CFL=$cfl" tools/run_batch.sh clainshock4 "RK4MUSCL4(EPD1);RK4MUSCL4(EPD2)"
done
SIM_SET="Time_CFL=0.8" tools/run_batch.sh clainshock4 "RK4MUSCL4(EPD1)"
