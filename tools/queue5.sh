#!/usr/bin/env bash
# Queue 5 (after queue 3): where does the 2-D shock CFL limit lie?
cd "$(dirname "$0")/.."
while pgrep -f "queue3.sh" > /dev/null; do sleep 60; done
M4="RK2MUSCL2(VK);RK2MUSCL2(EPD1);RK2MUSCL2(EPD2);RK2Upwind1"
for cfl in 2.5 3.0 4.0; do
  SIM_SET="Grid_Ns=(400,400); Time_CFL=$cfl" tools/run_batch.sh Euler2D "$M4"
done
