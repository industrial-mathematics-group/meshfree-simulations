#!/usr/bin/env bash
# Queue 4 (after queue 3): 2-D advection convergence with a narrower Gaussian (width 1: periodic-smooth to 1e-11;
# width 2 leaves a derivative kink of 5e-3 at the periodic boundary)
cd "$(dirname "$0")/.."
while pgrep -f "queue2.sh|queue3.sh" > /dev/null; do sleep 60; done
export SIM_PROGRESS=300
SIM_T=Float64 SIM_SET="IC_width=1.0" tools/run_batch.sh ConvergenceAdvection "RK4MUSCL2;RK4MUSCL3;RK4MUSCL4;RK4MUSCL5"
