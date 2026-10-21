#!/bin/bash
# Prints what the GPU node and the container see (driver, GPUs, CPU).
hostname; nvidia-smi; nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv; nproc; ls /usr/lib/x86_64-linux-gnu/libcuda* /.singularity.d/libs 2>&1 | head
