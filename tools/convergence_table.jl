# Convergence table from saved runs: reads every *.jld2 under <dir>/data (the save path layout), takes the
# final L1 error (stat :l1error) and the grid size of each run, and prints per experiment and method the errors
# and experimental orders EOC = log(e_{i-1}/e_i) / log(N_i/N_{i-1}). Experiments are told apart by the run
# parameters (dimension, PDE), not by file names.
#
#     julia --project=env_headless tools/convergence_table.jl [save dir, default .] [table.md] [table.csv]
using HyperCloud, StaticArrays, Printf
const JLD2 = HyperCloud.PDEStudioCore.JLD2

dir = joinpath(get(ARGS, 1, "."), "data")
runs = Dict{Tuple{String, String}, Dict{Int, Tuple{Float64, String}}}()   # (experiment, method + parameters) => N => (error, T)
for f in sort(filter(endswith(".jld2"), readdir(dir; join = true)))
    file = JLD2.load(f)
    p, st = file["params"], file["raw"].stats
    haskey(st, :l1error) || continue
    exp = "$(p[:PDE_D])-D $(p[:PDE_name])" * (haskey(p, :tmax) ? " (t = $(p[:tmax]))" : "")
    method = "$(p[:Time_stepper])$(p[:Scheme_name])$(p[:Scheme_order])" *
             (get(p, :MOOD_criterion, :none) == :none ? "" : "($(p[:MOOD_strategy]))") *
             (get(p, :Limiter_name, :none) == :none ? "" : "($(p[:Limiter_name]))")
    N = p[:Grid_Ns][1]
    err = last(st[:l1error])[1]
    # the point-cloud parameters that vary between sweeps (jitter, interpolation range, flux) are part of the key
    method *= " jitter $(p[:Grid_randomness_factor][1]) R $(p[:Weight_range])" *
              (get(p, :Flux_name, :Rusanov) == :Rusanov ? "" : " $(p[:Flux_name])")
    d = get!(runs, (exp, method), Dict{Int, Tuple{Float64, String}}())
    haskey(d, N) && @warn "duplicate run" exp method N
    d[N] = (err, string(p[:T]))
end

md = IOBuffer(); csv = IOBuffer()
println(csv, "experiment,method,N,T,l1error,eoc")
for (exp, method) in sort(collect(keys(runs)))
    d = runs[(exp, method)]; Ns = sort(collect(keys(d)))
    println(md, "\n### $exp, $method ($(d[Ns[1]][2]), $(length(Ns)) runs)\n")
    println(md, "| N | L1 error | EOC |\n|---|---|---|")
    for (i, N) in enumerate(Ns)
        eoc = i == 1 ? NaN : log(d[Ns[i-1]][1] / d[N][1]) / log(N / Ns[i-1])
        println(md, "| $N | ", @sprintf("%.3e", d[N][1]), " | ", i == 1 ? "–" : @sprintf("%.2f", eoc), " |")
        println(csv, exp, ",", method, ",", N, ",", d[N][2], ",", d[N][1], ",", i == 1 ? "" : eoc)
    end
end
# overview: EOC of the last three refinements per (experiment, method)
summary = IOBuffer()
println(summary, "| Experiment | Method | runs (N) | error at N_max | mean EOC (last 3) |\n|---|---|---|---|---|")
for (exp, method) in sort(collect(keys(runs)))
    d = runs[(exp, method)]; Ns = sort(collect(keys(d)))
    eocs = [log(d[Ns[i-1]][1] / d[Ns[i]][1]) / log(Ns[i] / Ns[i-1]) for i in max(2, length(Ns) - 2):length(Ns)]
    println(summary, "| $exp | $method | $(length(Ns)) ($(Ns[1])–$(Ns[end])) | ", @sprintf("%.2e", d[Ns[end]][1]), " | ",
            @sprintf("%.2f", sum(eocs) / length(eocs)), " |")
end
out = "# Convergence table ($(dir), $(sum(length, values(runs))) runs)\n\n" * String(take!(summary)) * String(take!(md))
print(out)
length(ARGS) >= 2 && write(ARGS[2], out)
length(ARGS) >= 3 && write(ARGS[3], String(take!(csv)))
