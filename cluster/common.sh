# Sourced by the job scripts.
exec 2>&1   # the gateway log only captures stdout
# Compiled package images are CPU specific: a depot per CPU type in front of the shared /depot keeps
# EPYC and Xeon nodes from overwriting each other's images (packages and artifacts stay in /depot).
cpu=$(julia --startup-file=no -e 'print(Sys.CPU_NAME)')
mkdir -p "/depot/cpu/$cpu"
export JULIA_DEPOT_PATH="/depot/cpu/$cpu:/depot"

# Method keys of the 2026-10-26 setup (experiments/methods_review.jl): k<K>.<variant>.<ar|ru>.mls<m>, e.g.
# k4.epd2.ar.mls2 -> "RK4SSPMUSCL4(EPD2,AR,MLS2)"; variant: none, epd1, epd2, vk, mm. Other keys pass through
# the per-script tables.
review_method() {
    local key=$1 K v f m st
    [[ $key =~ ^k([1-5])\.(none|epd1|epd2|vk|mm)\.(ar|ru)\.mls([0-9])$ ]] || return 1
    K=${BASH_REMATCH[1]}; v=${BASH_REMATCH[2]}; f=${BASH_REMATCH[3]}; m=${BASH_REMATCH[4]}
    case $v in epd1) v=EPD1 ;; epd2) v=EPD2 ;; vk) v=VK ;; mm) v=minmod ;; esac
    case $K in 1|2) st=RK2 ;; 3) st=RK3 ;; *) st=RK4SSP ;; esac
    echo "${st}MUSCL${K}(${v},${f^^},MLS${m})"
}
