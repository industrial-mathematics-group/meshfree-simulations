#!/usr/bin/env python3
"""EOC table from run_sim.jl summary lines: tools/eoc.py LOG [LOG...]  (L1 error of u1 at the final time vs N)"""
import re, sys, math, collections
runs = collections.defaultdict(dict)
for f in sys.argv[1:]:
    for line in open(f):
        m = re.match(r"(OK|FAIL)\s+(\S+) N=(\d+)(?:x\d+)*\s.*L1err\(u1\)=(\S+)", line)
        if m: runs[m.group(2)][int(m.group(3))] = float(m.group(4))
for meth, d in sorted(runs.items()):
    print(f"--- {meth}")
    Ns = sorted(d); prev = None
    for N in Ns:
        e = d[N]
        eoc = "" if prev is None else f"{-math.log(e/prev[1])/math.log(N/prev[0]):6.2f}"
        print(f"  N={N:6d}  L1={e:.3e}  EOC={eoc}")
        prev = (N, e)
