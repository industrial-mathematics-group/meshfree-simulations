#!/usr/bin/env bash
# Queue 2 (after queue 1): CFL study on the shock problems (GPU, Float32), then the 2-D convergence sweeps (GPU, Float64)
cd "$(dirname "$0")/.."
while pgrep -f queue1.sh > /dev/null; do sleep 60; done
export SIM_PROGRESS=300
M4="RK2MUSCL2(VK);RK2MUSCL2(EPD1);RK2MUSCL2(EPD2);RK2Upwind1"
for cfl in 0.2 0.4 0.6 0.8 1.0; do
  SIM_SET="Grid_Ns=(400,400); Time_CFL=$cfl" tools/run_batch.sh Euler2D "$M4"
done
for cfl in 0.2 0.4 0.6 0.8 1.0; do
  SIM_SET="Grid_Ns=(960,300); Time_CFL=$cfl" tools/run_batch.sh DMR "$M4"
done
for cfl in 0.2 0.4 0.6 0.8; do
  SIM_SET="Grid_Ns=(600,200); Time_CFL=$cfl" tools/run_batch.sh Mach3Step "RK2MUSCL2(VK);RK2MUSCL2(EPD1)"
done
SIM_T=Float64 tools/run_batch.sh ConvergenceAdvection "RK4MUSCL2;RK4MUSCL3;RK4MUSCL4;RK4MUSCL5"
SIM_T=Float64 tools/run_batch.sh ConvergenceBurgers "RK4MUSCL2;RK4MUSCL3;RK4MUSCL4;RK4MUSCL5"
