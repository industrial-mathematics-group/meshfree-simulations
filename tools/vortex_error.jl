# L1 error of the isentropic-vortex runs (experiments/Euler/IsentropicVortex.jl) against the initial condition at
# the final time (one period, so the exact solution is the IC), per method and N, with the order; and the mass
# conservation defect. Reads every *.jld2 under <dir>/data.
#     julia --project=env_headless tools/vortex_error.jl data/vortex [table.md]
using HyperCloud, StaticArrays, Printf
const JLD2 = HyperCloud.PDEStudioCore.JLD2
Base.include(Main, joinpath(@__DIR__, "..", "experiments", "Euler", "IsentropicVortex.jl"))

dir = joinpath(get(ARGS, 1, "data/vortex"), "data")
runs = Dict{Tuple{String, String}, Dict{Int, NamedTuple}}()
for f in sort(filter(endswith(".jld2"), readdir(dir; join = true)))
    file = JLD2.load(f); p, d = file["params"], file["raw"]
    get(p, :IC_name, :none) == :my_ic || continue
    label = "$(p[:Time_stepper])MUSCL$(p[:Scheme_order])" * (get(p, :MOOD_criterion, :none) == :none ? "(none" : "($(p[:MOOD_strategy])") *
            ",$(p[:Flux_name] == :AlignedRusanov ? "AR" : "RU"),MLS$(get(p, :Scheme_MLS_order, "-")))"
    x, u = d.x[end], d.u[end]
    vol = length(d.params) > 0 && haskey(file, "volumes") ? file["volumes"] : nothing
    n = length(u)
    # L1 error of the density against the IC (mean over particles, times the domain area 100 = integral L1 as elsewhere)
    err = sum(abs(u[i][1] - vortex_state(x[i][1], x[i][2])[1]) for i in 1:n) / n * 100
    m = [v[1] for v in d.stats[:mass]]
    runs[(string(p[:T]), label)] = get(runs, (string(p[:T]), label), Dict{Int, NamedTuple}())
    runs[(string(p[:T]), label)][p[:Grid_Ns][1]] = (err = err, defect = (m[end] - m[1]) / m[1], minrho = minimum(v -> v[1], u))
end
io = IOBuffer()
println(io, "| T | Method | N | L1(ρ) vs IC at t = 10 | order | mass defect | min ρ |\n|---|---|---|---|---|---|---|")
for k in sort(collect(keys(runs)))
    d = runs[k]; Ns = sort(collect(keys(d)))
    for (i, N) in enumerate(Ns)
        o = i == 1 ? "–" : @sprintf("%.2f", log(d[Ns[i-1]].err / d[N].err) / log(N / Ns[i-1]))
        println(io, "| $(k[1]) | $(k[2]) | $N | ", @sprintf("%.3e", d[N].err), " | $o | ", @sprintf("%+.2e", d[N].defect), " | ", @sprintf("%.4f", d[N].minrho), " |")
    end
end
out = String(take!(io)); print(out); length(ARGS) >= 2 && write(ARGS[2], out)
