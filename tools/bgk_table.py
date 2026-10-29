#!/usr/bin/env python3
"""Tables of the BGK-on-moving-grids study from the RESULT lines of the job logs (cluster/run_bgk.sh wave/sodc/vortex).

    python3 tools/bgk_table.py logs.txt > table.md

Tags: bgk<case>_<method>_<mover>_<vol>_<on|off>_tau<tau>[_p<on|off>]_N<n>. Per (case, method, mover, volumes,
adaptivity, tau, flux patch) and N: the L1 density error against the exact solution (wave, vortex), the relative
change of the interior content per moment (mass, momentum..., energy), the largest remap residual, the particle
count and min ρ. The order is the mean EOC of the last three refinements."""
import re, sys, math, collections

runs = collections.defaultdict(dict)
for line in open(sys.argv[1]):
    m = re.match(r"(?:\d+ )?RESULT tag=(\S+) (.*)", line.strip())
    if not m: continue
    tag, rest = m.groups()
    t = re.match(r"bgk(wave|sod|vortex)_(\w+?)_(none|lagrangian)_(\w+?)_(on|off)_tau([^_]+)(?:_p(on|off))?_N(\d+)$", tag)
    if not t: continue
    key = t.groups()[:7]; N = int(t.group(8))
    d = runs[key].setdefault(N, {})
    for k, v in re.findall(r"(\w+)=(\S+)", rest): d[k] = v

def comps(x): return [float(v) for v in x.split(",")]
def eoc(Ns, e):
    r = [math.log(e[Ns[i-1]] / e[Ns[i]]) / math.log(Ns[i] / Ns[i-1]) for i in range(1, len(Ns)) if e[Ns[i]] > 0 and e[Ns[i-1]] > 0]
    return sum(r[-3:]) / len(r[-3:]) if r else float("nan")
fmt = lambda x: "–" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.2f}"

out = ["# BGK on moving grids: accuracy and conservation", ""]
for case in ("wave", "sod", "vortex"):
    keys = sorted(k for k in runs if k[0] == case)
    if not keys: continue
    err = {"wave": "wave_L1_rho", "vortex": "vortex_L1_rho"}.get(case)
    out += [f"\n## {case}\n", "| method | mover | volumes | adapt | τ | ALE flux patch | N | " + ("L1(ρ) at N_max | order | " if err else "") +
            "interior change ρ / ρu / E at N_max | max remap | particles (+ins/−rem) | min ρ | status |",
            "|---|---|---|---|---|---|---|" + ("---|---|" if err else "") + "---|---|---|---|---|"]
    for k in keys:
        _, meth, mover, vol, ad, tau, patch = k
        d = runs[k]; Ns = sorted(n for n in d if "status" in d[n])
        if not Ns: continue
        n = Ns[-1]; dn = d[n]
        bad = [m for m in Ns if d[m].get("status") != "ok"]
        row = f"| {meth} | {mover} | {vol} | {ad} | {tau} | {patch or '–'} | {Ns[0]}–{n} ({len(Ns)}) | "
        if err:
            e = {m: float(d[m][err]) for m in Ns if err in d[m]}
            eN = sorted(e)
            row += (f"{e[eN[-1]]:.3e} | {fmt(eoc(eN, e))} | " if eN else "– | – | ")
        ch = comps(dn["rel_interior_change"]) if "rel_interior_change" in dn else []
        row += (" / ".join(f"{c:+.1e}" for c in (ch[:1] + ch[-1:] if len(ch) > 2 and case != "sod" else ch)) if ch else "–") + " | "
        rm = max((abs(c) for m in Ns if "max_abs_remap" in d[m] for c in comps(d[m]["max_abs_remap"])), default=0.0)
        row += f"{rm:.1e} | {dn.get('N', '?')} (+{dn.get('inserted', '?')}/−{dn.get('removed', '?')}) | {float(dn.get('min_u1', 'nan')):.3f} | "
        row += ("ok" if not bad else "failed at N=" + ",".join(map(str, bad))) + " |"
        out.append(row)
print("\n".join(out))
