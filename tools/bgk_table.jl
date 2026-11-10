# Tables of the BGK-on-moving-grids study from the RESULT lines of the job logs (cluster/run_bgk.sh wave/sodc/vortex),
# the Julia version of tools/bgk_table.py (same output; for containers without Python).
#
#     julia tools/bgk_table.jl logs.txt > table.md
#
# Tags: bgk<case>_<method>_<mover>_<vol>_<on|off>_tau<tau>[_p<on|off>]_N<n>. Per (case, method, mover, volumes, adaptivity,
# tau, flux patch) and N: the L1 density error against the exact solution (wave, vortex), the relative change of the
# interior content per moment (mass, momentum..., energy), the largest remap residual, the particle count and min ρ.
# The order is the mean EOC of the last three refinements. If a tag appears several times (a requeued job), the last
# line wins.
using Printf

const TAG = r"^bgk(wave|sod|vortex)_(\w+?)_(none|lagrangian)_(\w+?)_(on|off)_tau([^_]+)(?:_p(on|off))?_N(\d+)$"
runs = Dict{NTuple{7, String}, Dict{Int, Dict{String, String}}}()
for line in eachline(ARGS[1])
    m = match(r"^(?:\d+ )?RESULT tag=(\S+) (.*)$", strip(line))
    m === nothing && continue
    t = match(TAG, m.captures[1])
    t === nothing && continue
    key = ntuple(i -> something(t.captures[i], ""), 7)
    d = get!(get!(runs, key, Dict{Int, Dict{String, String}}()), parse(Int, t.captures[8]), Dict{String, String}())
    for kv in eachmatch(r"(\w+)=(\S+)", m.captures[2])
        d[kv.captures[1]] = kv.captures[2]
    end
end

comps(x) = parse.(Float64, split(x, ","))
function eoc(Ns, e)
    r = [log(e[Ns[i-1]] / e[Ns[i]]) / log(Ns[i] / Ns[i-1]) for i in 2:length(Ns) if e[Ns[i]] > 0 && e[Ns[i-1]] > 0]
    isempty(r) && return NaN
    last3 = r[max(1, end - 2):end]
    return sum(last3) / length(last3)
end
fmt(x) = isnan(x) ? "–" : @sprintf("%.2f", x)
sgn(c) = @sprintf("%+.1e", c)

out = ["# BGK on moving grids: accuracy and conservation", ""]
for case in ("wave", "sod", "vortex")
    keys_ = sort([k for k in keys(runs) if k[1] == case])
    isempty(keys_) && continue
    err = Dict("wave" => "wave_L1_rho", "vortex" => "vortex_L1_rho")
    e_key = get(err, case, nothing)
    push!(out, "\n## $case\n")
    push!(out, "| method | mover | volumes | adapt | τ | ALE flux patch | N | " * (e_key === nothing ? "" : "L1(ρ) at N_max | order | ") *
               "interior change ρ / ρu / E at N_max | max remap | particles (+ins/−rem) | min ρ | status |")
    push!(out, "|---|---|---|---|---|---|---|" * (e_key === nothing ? "" : "---|---|") * "---|---|---|---|---|")
    for k in keys_
        _, meth, mover, vol, ad, tau, patch = k
        d = runs[k]; Ns = sort([n for n in keys(d) if haskey(d[n], "status")])
        isempty(Ns) && continue
        n = Ns[end]; dn = d[n]
        bad = [m for m in Ns if get(d[m], "status", "") != "ok"]
        row = "| $meth | $mover | $vol | $ad | $tau | $(isempty(patch) ? "–" : patch) | $(Ns[1])–$n ($(length(Ns))) | "
        if e_key !== nothing
            e = Dict(m => parse(Float64, d[m][e_key]) for m in Ns if haskey(d[m], e_key))
            eN = sort(collect(keys(e)))
            row *= isempty(eN) ? "– | – | " : @sprintf("%.3e | %s | ", e[eN[end]], fmt(eoc(eN, e)))
        end
        ch = haskey(dn, "rel_interior_change") ? comps(dn["rel_interior_change"]) : Float64[]
        shown = length(ch) > 2 && case != "sod" ? [ch[1], ch[end]] : ch
        row *= (isempty(ch) ? "–" : join(sgn.(shown), " / ")) * " | "
        rm = maximum((abs(c) for m in Ns if haskey(d[m], "max_abs_remap") for c in comps(d[m]["max_abs_remap"])); init = 0.0)
        row *= @sprintf("%.1e | %s (+%s/−%s) | %.3f | ", rm, get(dn, "N", "?"), get(dn, "inserted", "?"), get(dn, "removed", "?"),
                        parse(Float64, get(dn, "min_u1", "NaN")))
        row *= (isempty(bad) ? "ok" : "failed at N=" * join(bad, ",")) * " |"
        push!(out, row)
    end
end
println(join(out, "\n"))
