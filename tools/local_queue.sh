#!/bin/bash
# Runs `cluster/run_method.sh`-style lines locally, one after the other (CPU, the precision of the line, 12 threads),
# saving into data/ like the cluster runs after collection. Each line: CONFIG KEY NX NY CFL T, KEY as in
# cluster/run_method.sh (m1, k<K>.<variant>.<ar|ru>.mls<m>). A line already marked done in the log is skipped, so the
# script is safe to rerun. Summary lines (OK/FAIL) of every run go to <log>.
#     setsid nohup bash tools/local_queue.sh <queue file> <log> &
queue=$1 log=$2
cd "$(dirname "$0")/.."
label() {   # same mapping as review_method in cluster/common.sh
    local key=$1 K v f m st
    [ "$key" = m1 ] && { echo RK2MUSCL1; return; }
    [[ $key =~ ^k([1-5])\.(none|epd1|epd2|vk|mm)\.(ar|ru)\.mls([0-9])$ ]] || return 1
    K=${BASH_REMATCH[1]}; v=${BASH_REMATCH[2]}; f=${BASH_REMATCH[3]}; m=${BASH_REMATCH[4]}
    case $v in epd1) v=EPD1 ;; epd2) v=EPD2 ;; vk) v=VK ;; mm) v=minmod ;; esac
    case $K in 1|2) st=RK2 ;; 3) st=RK3 ;; *) st=RK4SSP ;; esac
    echo "${st}MUSCL${K}(${v},${f^^},MLS${m})"
}
while read -r cfg key nx ny cfl T; do
    [[ -z $cfg || $cfg == \#* ]] && continue
    grep -qxF "# done $cfg $key $nx $ny $cfl $T" "$log" 2>/dev/null && continue
    method=$(label "$key") || { echo "unknown key $key" >> "$log"; continue; }
    echo "# start $cfg $key $nx $ny $cfl $T $(date +%F_%H:%M)" >> "$log"
    SIM_BACKEND=cpu SIM_T=$T SIM_SAVE=. SIM_METHODS="$method" SIM_SET="Grid_Ns=($nx,$ny); Time_CFL=$cfl" SIM_PROGRESS=600 \
        timeout 4h julia --project=env_headless -t 12 tools/run_sim.jl "$cfg" 2>&1 | grep -E '^ *(OK|FAIL)' | sed "s/^/$cfg $key $T: /" >> "$log"
    find data -maxdepth 1 -type f ! -perm 660 -exec chmod 660 {} +
    echo "# done $cfg $key $nx $ny $cfl $T" >> "$log"
done < "$queue"
echo "# queue finished $(date +%F_%H:%M)" >> "$log"
