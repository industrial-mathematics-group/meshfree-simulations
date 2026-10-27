#!/usr/bin/env python3
"""Tables of the particle-management study from the RESULT lines of the job logs (cluster/run_adapt.sh).

    python3 tools/adapt_table.py logs.txt > table.md

Tag format: <case>_K<k>_<vol>_<on|off>_g<maxgap>-<mingap>_s<spd>_e<excl>_N<n>. Per (case, K, vol, adaptivity,
settings): the relative mass change over the run (component 1, i.e. mass/density) per N with its order, the
remap's share (sum of the per-event mass residuals) and the particle count change."""
import re, sys, math, collections

runs = collections.defaultdict(dict)
for line in open(sys.argv[1]):
    m = re.match(r"RESULT tag=(\S+) (.*)", line.strip())
    if not m: continue
    tag, rest = m.groups()
    t = re.match(r"(\w+?)_K(\d)_(\w+?)_(on|off)_g([\d.]+)-([\d.]+)_s(\d+)_e(\d)_N(\d+)$", tag)
    if not t: continue
    key = t.groups()[:8]; N = int(t.group(9))
    d = runs[key].setdefault(N, {})
    for k, v in re.findall(r"(\w+)=(\S+)", rest): d[k] = v

def first(x): return float(x.split(",")[0])
out = []
summary = ["| case | K | volumes | adaptivity | settings (max/min gap, spd, excl) | N | status | mass change at N_max | order (last 3) | remap share (max) | particles N0→N_end (+ins/−rem) at N_max |",
           "|---|---|---|---|---|---|---|---|---|---|---|"]
for key in sorted(runs, key=lambda k: (k[0], int(k[1]), k[2], k[3], k[4:])):
    case, K, vol, ad, mg, ng, spd, ex = key
    d = runs[key]; Ns = sorted(n for n in d if "rel_mass_change" in d[n])
    if not Ns: continue
    dm = {n: first(d[n]["rel_mass_change"]) for n in Ns}
    eoc = [math.log(abs(dm[Ns[i-1]]) / abs(dm[Ns[i]])) / math.log(Ns[i] / Ns[i-1]) for i in range(1, len(Ns)) if dm[Ns[i]] != 0 and dm[Ns[i-1]] != 0]
    last = eoc[-3:]
    remap = max(abs(first(d[n].get("rel_remap_sum", "0"))) for n in Ns)
    n = Ns[-1]; st = d[n].get("status", "?")
    fails = [m for m in Ns if d[m].get("status", "ok") != "ok"]
    summary.append(f"| {case} | {K} | {vol} | {ad} | {mg}/{ng}, {spd}, {ex} | {Ns[0]}–{Ns[-1]} ({len(Ns)}) | {st}{' (fails at N=' + ','.join(map(str, fails)) + ')' if fails else ''} | "
                   f"{dm[n]:+.2e} | {sum(last)/len(last):.2f} | {remap:.1e} | {d[n].get('N', '?')} (+{d[n].get('inserted', '?')}/−{d[n].get('removed', '?')}) |" if last else
                   f"| {case} | {K} | {vol} | {ad} | {mg}/{ng}, {spd}, {ex} | {Ns[0]}–{Ns[-1]} ({len(Ns)}) | {st} | {dm[n]:+.2e} | – | {remap:.1e} | {d[n].get('N', '?')} |")
    out.append(f"\n### {case}, K = {K}, {vol} volumes, adaptivity {ad} ({mg}/{ng}, spd {spd}, excl {ex})\n")
    out.append("| N | status | rel. mass change (all components) | order | remap sum | max |Δm| per call | N0→N_end | +ins/−rem | min u₁ | runtime |")
    out.append("|---|---|---|---|---|---|---|---|---|---|")
    for i, m in enumerate(Ns):
        r = d[m]; o = "–" if i == 0 or dm[m] == 0 or dm[Ns[i-1]] == 0 else f"{math.log(abs(dm[Ns[i-1]]) / abs(dm[m])) / math.log(m / Ns[i-1]):.2f}"
        out.append(f"| {m} | {r.get('status','?')} | {r.get('rel_mass_change')} | {o} | {r.get('rel_remap_sum')} | {r.get('max_abs_remap')} | {r.get('N')} | +{r.get('inserted')}/−{r.get('removed')} | {r.get('min_u1')} | {r.get('runtime')} |")
print("# Particle-management study: mass conservation\n")
print("\n".join(summary)); print("\n".join(out))
