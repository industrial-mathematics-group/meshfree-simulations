#!/usr/bin/env python3
"""Side-by-side table of two convergence CSVs (from tools/convergence_table.jl), e.g. lattice vs jittered points:

    python3 tools/compare_convergence.py data/convergence/convergence_table.csv data/convergence_jitter02/convergence_table.csv

Prints, per experiment and method, the errors of both sets at every N and the ratio, plus the mean EOC of the last
three refinements of each set."""
import csv, math, sys
from collections import defaultdict

def load(path):
    d = defaultdict(dict)
    with open(path) as f:
        for r in csv.DictReader(f):
            d[(r["experiment"], r["method"])][int(r["N"])] = float(r["l1error"])
    return d

a, b = load(sys.argv[1]), load(sys.argv[2])
la, lb = sys.argv[3] if len(sys.argv) > 3 else "A", sys.argv[4] if len(sys.argv) > 4 else "B"
print(f"| Experiment | Method | N | {la} | {lb} | {lb}/{la} |\n|---|---|---|---|---|---|")
summary = []
for key in sorted(set(a) | set(b)):
    Ns = sorted(set(a.get(key, {})) | set(b.get(key, {})))
    for N in Ns:
        ea, eb = a.get(key, {}).get(N), b.get(key, {}).get(N)
        fmt = lambda e: f"{e:.3e}" if e is not None else "–"
        ratio = f"{eb / ea:.2f}" if ea and eb else "–"
        print(f"| {key[0]} | {key[1]} | {N} | {fmt(ea)} | {fmt(eb)} | {ratio} |")
    def eoc(d):
        n = sorted(k for k in Ns if k in d)
        if len(n) < 2: return float("nan")
        e = [math.log(d[n[i-1]] / d[n[i]]) / math.log(n[i] / n[i-1]) for i in range(max(1, len(n) - 3), len(n))]
        return sum(e) / len(e)
    summary.append((key, eoc(a.get(key, {})), eoc(b.get(key, {}))))
print(f"\n| Experiment | Method | mean EOC {la} (last 3) | mean EOC {lb} (last 3) |\n|---|---|---|---|")
for key, ea, eb in summary:
    print(f"| {key[0]} | {key[1]} | {ea:.2f} | {eb:.2f} |")
