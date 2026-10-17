#!/usr/bin/env bash
# Queue 3 (after queue 2 and cs4b): CFL limit beyond 1 on the shock problems, the skipped Mach3Step production run
cd "$(dirname "$0")/.."
while pgrep -f "queue2.sh|queue_cs4b.sh" > /dev/null; do sleep 60; done
export SIM_PROGRESS=300
M4="RK2MUSCL2(VK);RK2MUSCL2(EPD1);RK2MUSCL2(EPD2);RK2Upwind1"
for cfl in 1.2 1.5 2.0; do
  SIM_SET="Grid_Ns=(400,400); Time_CFL=$cfl" tools/run_batch.sh Euler2D "$M4"
  SIM_SET="Grid_Ns=(960,300); Time_CFL=$cfl" tools/run_batch.sh DMR "$M4"
done
tools/run_batch.sh Mach3Step "RK2Upwind1"
