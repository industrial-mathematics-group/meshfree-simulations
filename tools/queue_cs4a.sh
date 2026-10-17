#!/usr/bin/env bash
# clainshock4 study A: CFL x precision at 100x100 (RK4 MUSCL4 EPD1, default R=5.5)
cd "$(dirname "$0")/.."
export SIM_PROGRESS=300 SIM_SET="Grid_Ns=(100,100)"
for T in Float32 Float64; do
  for cfl in 0.01 0.05 0.1 0.2 0.4 0.8; do
    SIM_T=$T SIM_MSET="" SIM_SET="Grid_Ns=(100,100); Time_CFL=$cfl" tools/run_batch.sh clainshock4 "RK4MUSCL4(EPD1)"
  done
done
