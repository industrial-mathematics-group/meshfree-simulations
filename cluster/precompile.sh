#!/bin/bash
# Precompiles both environments for the CPU type of the node it runs on (optional: the first job on a
# node type compiles anyway; this just does it up front):  slurm submit cluster/precompile.sh cpus=8 mem=16 time=30
set -eo pipefail
source "$(dirname "$0")/common.sh"
for env in env_headless env_cuda; do
    echo "== $env ($cpu)"
    julia --project=$env -e 'using Pkg; Pkg.precompile(); using HyperCloud; println("ok")'
done
