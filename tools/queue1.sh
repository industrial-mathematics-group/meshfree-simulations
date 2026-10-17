#!/usr/bin/env bash
# Overnight queue 1: production runs at the configured resolutions (GPU, Float32)
cd "$(dirname "$0")/.."
export SIM_PROGRESS=120
tools/run_batch.sh Euler2D "RK2MUSCL2(EPD1);RK2MUSCL2(EPD2);RK2Upwind1"
tools/run_batch.sh DMR "RK2MUSCL2(VK);RK2MUSCL2(EPD1);RK2MUSCL2(EPD2);RK2Upwind1"
tools/run_batch.sh Mach3Step "RK2MUSCL2(VK);RK2MUSCL2(EPD1);RK2MUSCL2(EPD2);RK2Upwind1"
