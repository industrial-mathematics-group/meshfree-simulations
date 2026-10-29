# Summary tables of the Euler shock runs from saved results: per experiment and grid, one row per method with the
# final time reached, the minimum pressure and density at the end, the mass change (m(t_end)/m(0) − 1) and its
# difference to a reference method, steps and runtime. Reads every *.jld2 under <dir>/data.
#     julia --project=env_headless tools/euler_table.jl [save dir, default .] [table.md] [reference label]
using HyperCloud, StaticArrays, Printf
const JLD2 = HyperCloud.PDEStudioCore.JLD2

dir = joinpath(get(ARGS, 1, "."), "data")
refname = get(ARGS, 3, "RK2MUSCL2(EPD1,AR,MLS2)")
label(p) = "$(p[:Time_stepper])$(p[:Scheme_name])$(p[:Scheme_order])" *
    (get(p, :Limiter_name, :none) != :none ? "($(p[:Limiter_name])" : get(p, :MOOD_criterion, :none) != :none ? "($(p[:MOOD_strategy])" : "(none") *
    "," * (get(p, :Flux_name, :Rusanov) == :AlignedRusanov ? "AR" : "RU") * ",MLS$(get(p, :Scheme_MLS_order, "-")))"
rows = Dict{Tuple{String, String}, Dict{String, NamedTuple}}()   # (experiment, grid) => label => row
for f in sort(filter(endswith(".jld2"), readdir(dir; join = true)))
    file = JLD2.load(f); p, d = file["params"], file["raw"]
    p[:PDE_name] == :euler || continue
    exp = string(get(p, :Grid_domain, :rect)) * (haskey(p, :IC_states) ? "/riemann" : "") * (get(p, :Grid_periodic, false) ? "/periodic" : "") *
          (haskey(p, :Grid_bc_map) ? "" : "/bc=$(get(p, :Grid_bc, "?"))")
    exp = haskey(p, :Grid_bc_map) ? exp : exp * "-outflow"
    grid = join(p[:Grid_Ns], "x") * " CFL $(p[:Time_CFL]) $(p[:T])"
    u = d.u[end]; γ = p[:PDE_gamma]
    pr(s) = (γ - 1) * (s[end] - sum(abs2, s[2:end-1]) / (2s[1]))
    m = [x[1] for x in d.stats[:mass]]
    st = d.stats
    get!(rows, (exp, grid), Dict{String, NamedTuple}())[label(p)] = (t = d.t[end], tmax = p[:tmax], minp = minimum(pr, u), minrho = minimum(s -> s[1], u),
        mass = m[end] / m[1] - 1, steps = haskey(st, :k_step) ? Int(st[:k_step][1]) : -1, runtime = haskey(st, :Runtime) ? st[:Runtime][1] : NaN)
end
io = IOBuffer()
for key in sort(collect(keys(rows)))
    r = rows[key]; ref = get(r, refname, nothing)
    println(io, "\n### $(key[1]), $(key[2])  (Δmass relative to the final mass of $refname)\n")
    println(io, "| Method | t_end | min ρ | min p | mass change | Δmass | steps | runtime |\n|---|---|---|---|---|---|---|---|")
    for lab in sort(collect(keys(r)))
        x = r[lab]; ok = x.t >= x.tmax * (1 - 1e-4) && x.minp > 0 && x.minrho > 0
        dm = ref === nothing ? "–" : @sprintf("%+.1f %%", 100 * (x.mass - ref.mass) / (1 + ref.mass))
        println(io, "| $lab | ", ok ? @sprintf("%.3g", x.t) : @sprintf("**%.3g (failed)**", x.t), " | ", @sprintf("%.3e", x.minrho), " | ", @sprintf("%.3e", x.minp),
                " | ", @sprintf("%+.4f", x.mass), " | $dm | $(x.steps) | ", @sprintf("%.0f s", x.runtime), " |")
    end
end
out = "# Euler shock problems ($(dir))\n" * String(take!(io)); print(out); length(ARGS) >= 2 && write(ARGS[2], out)
