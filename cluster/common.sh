# Sourced by the job scripts.
exec 2>&1   # the gateway log only captures stdout
# Compiled package images are CPU specific: a depot per CPU type in front of the shared /depot keeps
# EPYC and Xeon nodes from overwriting each other's images (packages and artifacts stay in /depot).
cpu=$(julia --startup-file=no -e 'print(Sys.CPU_NAME)')
mkdir -p "/depot/cpu/$cpu"
export JULIA_DEPOT_PATH="/depot/cpu/$cpu:/depot"
