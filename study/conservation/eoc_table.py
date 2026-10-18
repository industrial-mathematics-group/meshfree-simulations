#!/usr/bin/env python3
"""EOC table of a study TSV: eoc_table.py FILE KEYCOLS... --n N --err COL [--extra COL]"""
import csv, sys, math, collections, argparse
ap = argparse.ArgumentParser(); ap.add_argument("file"); ap.add_argument("keys", nargs="+")
ap.add_argument("--n", default="N"); ap.add_argument("--err", default="L1"); ap.add_argument("--extra", default=None)
a = ap.parse_args()
rows = list(csv.DictReader(open(a.file), delimiter="\t"))
t = collections.OrderedDict()
for r in rows: t.setdefault(tuple(r[k] for k in a.keys), {})[int(r[a.n])] = r
for k, v in t.items():
    Ns = sorted(v); e = [float(v[N][a.err]) for N in Ns]
    eoc = ["  -  "] + [f"{-math.log(e[i]/e[i-1])/math.log(Ns[i]/Ns[i-1]):5.2f}" if e[i] > 0 and e[i-1] > 0 else "  ?  " for i in range(1, len(Ns))]
    s = " ".join(f"{x:.2e}({o})" for x, o in zip(e, eoc))
    ex = "" if a.extra is None else "  " + a.extra + ": " + " ".join(f"{abs(float(v[N][a.extra])):.0e}" for N in Ns)
    print(" ".join(f"{x:>8}" for x in k), "|", s, ex)
