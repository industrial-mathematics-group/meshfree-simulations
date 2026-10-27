# Mass-conservation tables from saved runs: reads every *.jld2 under <dir>/data and, per run, the signed relative
# mass error at the end time (relative_mass[end] - 1 of the hyperbolic stat preset: numerical / reference mass)
# and its time integral (spacetime_relative_mass - 1). Groups by experiment (dimension, PDE), method, jitter,
# range, alpha, MLS order and N; prints one markdown table per (experiment, parameter set) with the error over N
# and its convergence order, plus a CSV of all runs.
#
#     julia --project=env_headless tools/mass_table.jl data/mass [tables.md] [runs.csv]
using HyperCloud, StaticArrays, Printf
const JLD2 = HyperCloud.PDEStudioCore.JLD2

dir = joinpath(get(ARGS, 1, "data/mass"), "data")
struct Key; exp::String; method::String; jitter::Float64; range::Float64; alpha::Float64; mls::Int; end
Base.isless(a::Key, b::Key) = (a.exp, a.method, a.jitter, a.range, a.alpha, a.mls) < (b.exp, b.method, b.jitter, b.range, b.alpha, b.mls)
runs = Dict{Key, Dict{Int, NamedTuple}}()
for f in sort(filter(endswith(".jld2"), readdir(dir; join = true)))
    file = JLD2.load(f); p, d = file["params"], file["raw"]; st = d.stats
    haskey(st, :relative_mass) || continue
    exp = "$(p[:PDE_D])-D $(p[:PDE_name])" * (get(p, :Grid_periodic, false) ? " periodic" : "") * " (t = $(p[:tmax]))"
    method = "$(p[:Time_stepper])$(p[:Scheme_name])$(p[:Scheme_order])" *
             (get(p, :MOOD_criterion, :none) == :none ? "" : "($(p[:MOOD_strategy]))") *
             (get(p, :Limiter_name, :none) == :none ? "" : "($(p[:Limiter_name]))")
    k = Key(exp, method, Float64(p[:Grid_randomness_factor][1]), Float64(p[:Weight_range]), Float64(p[:Weight_alpha]),
            Int(get(p, :Scheme_MLS_order, -1)))
    rm = st[:relative_mass]
    rel_end = last(rm)[1] - 1
    rel_max = maximum(x -> abs(x[1] - 1), rm)
    rel_st = haskey(st, :spacetime_relative_mass) ? st[:spacetime_relative_mass][1] - 1 : NaN
    l1 = haskey(st, :l1error) ? last(st[:l1error])[1] : NaN
    # Conservation defect from the discrete mass Σ V_i u_i itself (sampling-free): the change of the mass minus
    # the exact boundary flux (zero for periodic domains; (f(u_L) - f(u_R)) t for a Riemann problem with outflow
    # boundaries whose shock stays inside, f = u²/2 for Burgers, u·a for linear advection), relative to the mass.
    # Component 1 (mass; for Euler also momentum and energy below) of Σ V_i U_i over time
    mvec = st[:mass]; M = length(mvec[1])
    m = [x[1] for x in mvec]
    expected = zeros(M)
    if !get(p, :Grid_periodic, false) && get(p, :IC_name, :none) == :riemann
        UL, UR = collect(Float64, p[:IC_uL]), collect(Float64, p[:IC_uR])
        if p[:PDE_name] == :burgers
            expected[1] = (UL[1]^2 - UR[1]^2) / 2 * d.t[end]
        elseif p[:PDE_name] == :linear
            expected[1] = (UL[1] - UR[1]) * p[:PDE_velocities][1][1] * d.t[end]
        elseif p[:PDE_name] == :euler && M == 3      # 1-D Euler: F = (ρu, ρu² + p, u(E + p)), p = (γ-1)(E - ½ρu²)
            γ = p[:PDE_gamma]
            F(U) = (u = U[2] / U[1]; pr = (γ - 1) * (U[3] - U[2]^2 / (2U[1])); [U[2], U[2] * u + pr, u * (U[3] + pr)])
            expected .= (F(UL) .- F(UR)) .* d.t[end]
        end
    end
    defect = (m[end] - m[1] - expected[1]) / abs(m[end])
    # per-component defects (Euler: mass, momentum, energy), relative to max(|m_c(t_end)|, |m_c(0)|, 1e-300)
    defects = [(mvec[end][c] - mvec[1][c] - expected[c]) / max(abs(mvec[end][c]), abs(mvec[1][c]), 1e-300) for c in 1:M]
    get!(runs, k, Dict{Int, NamedTuple}())[p[:Grid_Ns][1]] = (rel_end = rel_end, rel_max = rel_max, rel_st = rel_st, l1 = l1, defect = defect, defects = defects, T = string(p[:T]))
end

md = IOBuffer(); csv = IOBuffer()
println(csv, "experiment,method,jitter,range,alpha,mls,N,T,conservation_defect,rel_mass_err_end,rel_mass_err_max,rel_mass_err_spacetime,l1error")
summary = IOBuffer()
println(summary, "| Experiment | Method | jitter | R | α | MLS | N range | conservation defect at N_max | order of the defect (last 3) | rel. mass error vs sampled reference (end / max) |")
println(summary, "|---|---|---|---|---|---|---|---|---|---|")
for k in sort(collect(keys(runs)))
    d = runs[k]; Ns = sort(collect(keys(d)))
    println(md, "\n### $(k.exp), $(k.method), jitter $(k.jitter), R = $(k.range), α = $(k.alpha), MLS $(k.mls) ($(d[Ns[1]].T))\n")
    println(md, "| N | conservation defect | order | rel. mass error (end) | max over t | space-time | L1 error |\n|---|---|---|---|---|---|---|")
    orders = Float64[]
    for (i, N) in enumerate(Ns)
        r = d[N]
        o = i == 1 ? NaN : log(abs(d[Ns[i-1]].defect) / abs(r.defect)) / log(N / Ns[i-1])
        i > 1 && push!(orders, o)
        comp = length(r.defects) > 1 ? " (components: " * join((@sprintf("%+.2e", x) for x in r.defects), ", ") * ")" : ""
        println(md, "| $N | ", @sprintf("%+.3e", r.defect), comp, " | ", i == 1 ? "–" : @sprintf("%.2f", o), " | ", @sprintf("%+.3e", r.rel_end), " | ",
                @sprintf("%.3e", r.rel_max), " | ", @sprintf("%+.3e", r.rel_st), " | ", @sprintf("%.3e", r.l1), " |")
        println(csv, join((k.exp, k.method, k.jitter, k.range, k.alpha, k.mls, N, r.T, r.defect, r.rel_end, r.rel_max, r.rel_st, r.l1), ","))
    end
    last3 = orders[max(1, end - 2):end]
    rN = d[Ns[end]]
    println(summary, "| $(k.exp) | $(k.method) | $(k.jitter) | $(k.range) | $(k.alpha) | $(k.mls) | $(Ns[1])–$(Ns[end]) ($(length(Ns))) | ",
            @sprintf("%+.2e", rN.defect), " | ", isempty(last3) ? "–" : @sprintf("%.2f", sum(last3) / length(last3)), " | ",
            @sprintf("%+.2e / %.2e", rN.rel_end, rN.rel_max), " |")
end
out = "# Mass conservation ($(dir), $(sum(length, values(runs))) runs)\n\nConservation defect = (m(t_end) − m(0) − exact boundary flux) / m(t_end) from the discrete mass Σ V_i u_i (sampling-free; the exact flux is 0 for periodic domains and (f(u_L) − f(u_R)) t_end for the outflow Riemann problem). The relative mass error of the hyperbolic stat preset (numerical / sampled reference − 1) is given for comparison; it contains the O(Δx) sampling error of a discontinuous reference. The order is that of |defect| over the last three refinements.\n\n" *
      String(take!(summary)) * String(take!(md))
print(out)
length(ARGS) >= 2 && write(ARGS[2], out)
length(ARGS) >= 3 && write(ARGS[3], String(take!(csv)))
