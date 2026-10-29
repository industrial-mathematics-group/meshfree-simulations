#!/usr/bin/env python3
"""Summary tables from the run summary lines (OK/FAIL) of the job logs, without the data files: the convergence
sweeps (cluster/run_mass.sh: L1 error per N and the mean EOC of the last three refinements), the Euler shock problems
(cluster/run_method.sh: final time, positivity, mass change, steps, runtime) and the vortex sweeps (cluster/run_sweep.sh).

    python3 tools/okline_tables.py <results.txt> <queue log>... > tables.md

results.txt: lines "<job id> <log line>" (cluster/queues/collect.sh); the queue logs map job ids to the queue lines,
whose arguments identify the run (experiment, order, variant, jitter, range, MLS degree)."""
import re, sys, math, collections

jobs = {}
for q in sys.argv[2:]:
    for l in open(q):
        m = re.match(r"(\d+) (.*?) :: (\S+) (.*)", l.strip())
        if m: jobs[m.group(1)] = (m.group(3).split("/")[-1], m.group(4).split(), m.group(2))

lines = collections.defaultdict(list)
for l in open(sys.argv[1]):
    m = re.match(r"(\d+) +((?:OK|FAIL) .*)", l.rstrip())
    if m and m.group(1) in jobs: lines[m.group(1)].append(m.group(2))

def field(s, k):
    m = re.search(rf"(?<![\w(]){re.escape(k)}=(\S+)", s)
    return m.group(1) if m else None
def num(x):
    try: return float(x)
    except (TypeError, ValueError): return float("nan")
def drift(s):
    m = re.search(r"mass drift=(\S+)", s); return m.group(1) if m else "–"
def N_of(s):
    m = re.search(r"N=(\d+)", s); return int(m.group(1)) if m else 0
def method_of(s):
    m = re.match(r"(?:OK|FAIL)\s+(\S+)", s); return m.group(1) if m else "?"

out = []
# Convergence: run_mass.sh NAME ORDER VARIANT JITTER RANGE ALPHA MLS NMIN NMAX
conv = collections.defaultdict(dict)
for j, (script, a, _) in jobs.items():
    if script != "run_mass.sh": continue
    key = (a[0], int(a[1]), a[2], a[3], a[4], a[6])
    for s in lines[j]:
        if s.startswith("OK"): conv[key][N_of(s)] = num(field(s, "L1err(u1)"))
        else: conv[key].setdefault(N_of(s), float("nan"))
if conv:
    out += ["## Convergence (CPU, Float64)", "",
            "| experiment | K | variant | jitter | R | MLS | N | L1 at N_max | EOC (mean of last 3) | failed |", "|---|---|---|---|---|---|---|---|---|---|"]
    for key in sorted(conv):
        d = conv[key]; Ns = sorted(n for n in d if d[n] == d[n] and d[n] > 0)
        bad = sorted(n for n in d if not (d[n] == d[n]))
        e = [math.log(d[Ns[i-1]] / d[Ns[i]]) / math.log(Ns[i] / Ns[i-1]) for i in range(1, len(Ns))]
        last = e[-3:]
        out.append(f"| {key[0]} | {key[1]} | {key[2]} | {key[3]} | {key[4]} | {key[5]} | {Ns[0] if Ns else '–'}–{Ns[-1] if Ns else '–'} ({len(Ns)}) | "
                   f"{d[Ns[-1]]:.2e} | " + (f"{sum(last)/len(last):.2f}" if last else "–") + f" | {','.join(map(str, bad)) or '–'} |")

# Euler shock problems: run_method.sh NAME KEY NX NY CFL BACKEND T
eul = []
for j, (script, a, res) in jobs.items():
    if script != "run_method.sh": continue
    for s in lines[j] or ["FAIL ? (no summary line)"]:
        eul.append((a[0], f"{a[2]}x{a[3]}", f"{a[5]} {a[6]}", f"{a[1]} ({method_of(s)})", s, j))
if eul:
    out += ["", "## Euler shock problems", "",
            "| experiment | grid | backend | method key (label) | status | t | NaN | min ρ | min p | mass change | steps | runtime | job |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for e, g, b, m, s, j in sorted(eul):
        out.append(f"| {e} | {g} | {b} | {m} | {s.split()[0]} | {field(s, 't') or '–'} | {field(s, 'NaN') or '–'} | {field(s, 'min(u1)') or '–'} | "
                   f"{field(s, 'min(p)') or '–'} | {drift(s)} | {field(s, 'steps') or '–'} | {field(s, 'runtime') or '–'} | {j} |")

# Vortex sweeps: run_sweep.sh NAME KEY NMIN NMAX BACKEND T
vor = []
for j, (script, a, res) in jobs.items():
    if script != "run_sweep.sh": continue
    for s in lines[j]:
        vor.append((a[1], f"{a[4]} {a[5]}", N_of(s), s, j))
if vor:
    out += ["", "## Isentropic vortex (end state; the error against the IC needs the data files, tools/vortex_error.jl)", "",
            "| method key | backend | N | status | NaN | min ρ | min p | mass change | steps | runtime |", "|---|---|---|---|---|---|---|---|---|---|"]
    for k, b, N, s, j in sorted(vor):
        out.append(f"| {k} | {b} | {N} | {s.split()[0]} | {field(s, 'NaN') or '–'} | {field(s, 'min(u1)') or '–'} | {field(s, 'min(p)') or '–'} | "
                   f"{drift(s)} | {field(s, 'steps') or '–'} | {field(s, 'runtime') or '–'} |")
print("\n".join(out))
